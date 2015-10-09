#define WEIGHTED 1
#include "ligra.h"


int K = 20; //dimensions of the latent vector

template <class vertex>
struct GD_F {
  double* latent_curr, *error;
  vertex* V;
  GD_F(double* _latent_curr, double* _error, vertex* _V) : 
    latent_curr(_latent_curr), error(_error), V(_V) {}
  inline bool update(uintE s, uintE d, intE edgeLen){ //update function applies PageRank equation
    double estimate = 0;
    for(int i = 0; i < K; i++){
      estimate += latent_curr[K*d+i]*latent_curr[K*s+i];
    }
    double err = edgeLen - estimate;

    for (int i = 0; i < K; i++){
      error[K*d + i] += latent_curr[K*s + i]*err;
    }
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen) { //atomic Update
    //no needed as we will always do pull based
    //writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct GD_Vertex_F {
  double step;
  double lambda;
  double* latent_curr;
  double* error;
  GD_Vertex_F(double* _latent_curr, double* _error, double _step, double _lambda) :
    latent_curr(_latent_curr), error(_error), 
    step(_step), lambda(_lambda){}
  inline bool operator () (uintE i) {
    for (int j = 0; j < K; j++){
      latent_curr[K*i + j] += step*(-lambda*latent_curr[K*i + j] + error[K*i + j]);
      error[K*i+j] = 0.0;
    }
    return 1;
  }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  
  //initialize latent vectors and errors
  srand(0);
  const intE n = GA.n;
  double* latent_curr = newA(double, K*n);
  double* error = newA(double, K*n);
  int numIter = 10;
  double step = 0.00000035;
  double lambda = 0.001;

  parallel_for(int i = 0; i < K*n; i++){
    latent_curr[i] = rand();
    error[i] = 0.0;
  }

  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);

  for (int iter = 0; iter < numIter; iter++){
    //edgemap to accumulate error for each node
    edgeMap(GA, Frontier, GD_F<vertex>(latent_curr,error,GA.V),GA.m/20);
    //vertexmap to update the latent vectors
    vertexMap(Frontier,GD_Vertex_F(latent_curr,error,step,lambda));
  }

  Frontier.del();

}
