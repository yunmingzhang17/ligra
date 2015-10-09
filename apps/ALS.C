#define WEIGHTED 1
#include "ligra.h"
#include "eigen_wrapper.hpp"
#include <stdlib.h>

#define COMPUTE_RMSE 1 //single threaded debug flag
//#define DEBUG 1
//#define PRECOMPUTE 1
int D = 20; //number of latent factors
double lambda = 0.065;
double minval = 1; //max allowed value in matrix                         
double maxval = 5; //min allowed value in matrix
int numUsers = 0;
int threshold = 0;

#ifdef COMPUTE_RMSE
double rmse;
#endif

struct vertex_data {
  vec pvec;

  vertex_data() {
    pvec = zeros(D);
  }
  void set_val(int index, float val){
    pvec[index] = val;
  }
  float get_val(int index){
    return pvec[index];
  }
};

std::vector<vertex_data> latent_factors_inmem;

#ifdef PRECOMPUTE
std::map<int, mat> largeOutDegreeMap;
void computeHighOutDegreeXTX(std::vector<int>* largeOutDegreeList){
  parallel_for (int i = 0; i < largeOutDegreeList->size(); i++){
    int vertexId = (*largeOutDegreeList)[i];
    largeOutDegreeMap[vertexId] = (latent_factors_inmem[vertexId].pvec*latent_factors_inmem[vertexId].pvec.transpose());
  }
}

#endif

/** compute a missing value based on ALS algorithm */
float als_predict(const vertex_data& user, 
		  const vertex_data& movie, 
		  const float rating, 
		  double & prediction, 
		  void * extra = NULL){


  prediction = dot_prod(user.pvec, movie.pvec);
  //truncate prediction to allowed values
  prediction = std::min((double)prediction, maxval);
  prediction = std::max((double)prediction, minval);
  //return the squared error
  float err = rating - prediction;
  assert(!std::isnan(err));
  return err*err; 

}


template <class vertex>
struct ALS_Vertex_F {
  vertex* V;
  //vector<float> * latent_factors_inmem;


  ALS_Vertex_F(vertex* _V) :
    V(_V){};
 
  inline bool operator() (uintE i){
#ifdef DEBUG
    cout << "vertex: " << i << endl;
    cout << "latent factor " << endl;
    for (int j = 0; j < D; j++){
      cout << latent_factors_inmem[i].pvec[j] << endl;
    }
    cout << endl;
#endif
    vertex_data & vdata = latent_factors_inmem[i];
    mat XtX = mat::Zero(D, D); 
    vec Xty = vec::Zero(D);
    uintE d = V[i].getInDegree();
    for (uintE j = 0; j < d; j++){
      uintE ngh = V[i].getInNeighbor(j);

      

      float observation = V[i].getInWeight(j);
      vertex_data & nbr_latent = latent_factors_inmem[ngh];
#ifdef DEBUG
      cout << "edge src: " << i << " dst: " << ngh << " rating: " << observation << endl;
      
#endif

      Xty += nbr_latent.pvec * observation;

#ifdef PRECOMPUTE    
      
      //std::map<int,mat>::iterator it = largeOutDegreeMap.find(ngh);
      //if(it != largeOutDegreeMap.end()){
      //if (largeOutDegreeMap.count(ngh)){
      if (V[ngh].getOutDegree() > threshold){
	  //element found;
	XtX.triangularView<Eigen::Upper>() += largeOutDegreeMap[ngh];
      } else {
	XtX.triangularView<Eigen::Upper>() += nbr_latent.pvec * nbr_latent.pvec.transpose();
      }
#else

      XtX.triangularView<Eigen::Upper>() += nbr_latent.pvec * nbr_latent.pvec.transpose();
      //nbr_latent.pvec.transpose();

#endif
      //XtX += nbr_latent.pvec * nbr_latent.pvec.transpose();

#ifdef COMPUTE_RMSE
      if (i < numUsers){//acumulate RMSE only for users
	double prediction;
	rmse += als_predict(vdata, nbr_latent, observation, prediction);
      }
#endif

    }    
    double regularization = lambda;
    //if (regnormal)
    //regularization *= vertex.num_edges();
    for(int i=0; i < D; i++) XtX(i,i) += regularization;
    //XtX = XtX.triangularView<Eigen::Upper>();
    vdata.pvec = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xty);

  }
};

#ifdef COMPUTE_RMSE

double training_rmse(int iteration, long numEdges ){  
  numEdges = numEdges/2; //because we are doubling the number of edges compare to Graphchi
  
  double ret = sqrt(rmse/numEdges); 
  cout << "Iteration: " << iteration << " rmse: " << rmse <<  " Training RMSE: " << ret << endl;
  return ret;
}

#endif

template<typename T>
void init_feature_vectors(uint size, T& latent_factors_inmem, bool randomize = true, double scale = 1.0){
  assert(size > 0);
  srand48(1);
  latent_factors_inmem.resize(size); // Initialize in-memory vertices.                 
  if (!randomize)
    return;
  for (int i=0; i < (int)size; i++){
    for (int j=0; j<D; j++)
      latent_factors_inmem[i].pvec[j] = scale * drand48();
  }
}
  
template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  numUsers = P.getOptionLongValue("-nusers",0);
  long n = GA.n;
  long numEdges = GA.m;


#ifdef COMPUTE_RMSE
  rmse = 0;
#endif

  cout << "num vertices: " << n << endl;
  cout << "num edges: " << numEdges << endl;


  
  init_feature_vectors<std::vector<vertex_data> >(n, latent_factors_inmem);


#ifdef PRECOMPUTE
  //do a sequential intiialization to avoid data races (should use a parallel version later
  startTime();
  threshold = 400;
  std::vector<int>* largeOutDegreeNodes = new std::vector<int>();
  for (long i = 0; i < n; i++) {
    if(GA.V[i].getOutDegree() > threshold){
      largeOutDegreeNodes->push_back(i);
      largeOutDegreeMap[i] = latent_factors_inmem[i].pvec*latent_factors_inmem[i].pvec.transpose(); 
    }
  }
  nextTime("Init time");
  startTime();
#endif

  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);

  for (int iter = 0; iter < 5; iter++){    
    vertexMap(Frontier, ALS_Vertex_F<vertex>(GA.V));
#ifdef PRECOMPUTE
    computeHighOutDegreeXTX(largeOutDegreeNodes);
#endif

#ifdef COMPUTE_RMSE
    training_rmse(iter, numEdges);
    rmse = 0;
#endif
  }
  Frontier.del();

}
