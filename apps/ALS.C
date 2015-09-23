#define WEIGHTED 1
#include "ligra.h"
#include "eigen_wrapper.hpp"
#include <stdlib.h>

#define COMPUTE_RMSE 1 //single threaded debug flag

int D = 20; //number of latent factors
double lambda = 0.065;
double minval = -1e100; //max allowed value in matrix                         
double maxval = 1e100; //min allowed value in matrix

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
  std::vector<vertex_data> latent_factors_inmem;

  ALS_Vertex_F(vertex* _V, std::vector<vertex_data> _latent_factors_inmem) :
    V(_V), latent_factors_inmem(_latent_factors_inmem){};
 
  inline bool operator() (uintE i){
#ifdef DEBUG
    cout << "vertex: " << i << endl;
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
      XtX.triangularView<Eigen::Upper>() += nbr_latent.pvec * nbr_latent.pvec.transpose();

#ifdef COMPUTE_RMSE
      double prediction;
      rmse += als_predict(vdata, nbr_latent, observation, prediction);
#endif

    }    
    double regularization = lambda;
    //if (regnormal)
    //regularization *= vertex.num_edges();
    for(int i=0; i < D; i++) XtX(i,i) += regularization;
    vdata.pvec = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xty);

  }
};


double training_rmse(int iteration, long numEdges ){  
  double ret = sqrt(rmse/numEdges); 
  cout << "Iteration: " << iteration << " Training RMSE: " << ret << endl;
  return ret;
}



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
  long n = GA.n;
  long numEdges = GA.m;
  rmse = 0;
  cout << "num vertices: " << n << endl;
  cout << "num edges: " << numEdges << endl;

  std::vector<vertex_data> latent_factors_inmem;
  init_feature_vectors<std::vector<vertex_data> >(n, latent_factors_inmem);

  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);
  vertexMap(Frontier, ALS_Vertex_F<vertex>(GA.V, latent_factors_inmem));
#ifdef COMPUTE_RMSE
  training_rmse(1, numEdges);
#endif
  rmse = 0;

  delete[] latent_factors_inmem;
  Frontier.del();
}
