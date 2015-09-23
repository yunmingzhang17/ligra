#define WEIGHTED 1
#include "ligra.h"
#include "eigen_wrapper.hpp"

int D = 20; //number of latent factors
double lambda = 0.065;


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

template <class vertex>
struct ALS_Vertex_F {
  vertex* V;
  //vector<float> * latent_factors_inmem;
  std::vector<vertex_data> latent_factors_inmem;

  ALS_Vertex_F(vertex* _V) :
    V(_V){};
 
  inline bool operator() (uintE i){
    cout << "vertex: " << i << endl;
    vertex_data & vdata = latent_factors_inmem[i];
    mat XtX = mat::Zero(D, D); 
    vec Xty = vec::Zero(D);


    uintE d = V[i].getInDegree();
    for (uintE j = 0; j < d; j++){
      uintE ngh = V[i].getInNeighbor(j);
      intE rating = V[i].getInWeight(j);
      cout << "edge src: " << i << " dst: " << ngh << " rating: " << rating << endl; 
    }
    
    double regularization = lambda;
    //if (regnormal)
    //regularization *= vertex.num_edges();
    for(int i=0; i < D; i++) XtX(i,i) += regularization;
    vdata.pvec = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xty);

  }
};
  
template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long n = GA.n;
  cout << "num vertices: " << n << endl;
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);
  vertexMap(Frontier, ALS_Vertex_F<vertex>(GA.V));
  Frontier.del();
}
