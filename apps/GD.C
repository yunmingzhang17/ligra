#define WEIGHTED 1
#include "ligra.h"
//#define COMPUTE_RMSE 1
//#define DEBUG2 1
//#define DEBUG 1

int K = 20; //dimensions of the latent vector

#ifdef COMPUTE_RMSE
double rmse;
#endif

template <class vertex>
struct GD_F {
  double* latent_curr, *error;
  vertex* V;
  GD_F(double* _latent_curr, double* _error, vertex* _V) : 
    latent_curr(_latent_curr), error(_error), V(_V) {}
  inline bool update(uintE s, uintE d, intE edgeLen){ //update function applies PageRank equation

#ifdef DEBUG
    cout << "edge src: " << s << " dst: " << d << " rating: " << edgeLen << endl;
#endif

    double estimate = 0;
    for(int i = 0; i < K; i++){
      estimate += latent_curr[K*d+i]*latent_curr[K*s+i];
    }
    double err = edgeLen - estimate;

    #ifdef COMPUTE_RMSE
    rmse += err*err;
    #endif

#ifdef DEBUG
    cout << "rating: " << edgeLen << endl;
    cout << "estimate: " << estimate << endl;
#endif

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

#ifdef DEBUG
    cout << "vertex: " << i << endl;
    cout << "latent_factor: " << endl;
    for (int j = 0; j < K; j++){
      cout << latent_curr[K*i+j] << endl;
    }
    cout << endl;
#endif

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
  int numIter = 5;
#ifndef DEBUG
  double step = 0.00000035;
#else
  double step = 0.0035;
#endif

  double lambda = 0.001;

#ifdef DEBUG2
  cout << "num vertices: " << n << endl;
  cout << "num edges: " << GA.m << endl;
#endif

  parallel_for(int i = 0; i < K*n; i++){
#ifndef COMPUTE_RMSE
    latent_curr[i] = ((double)rand()/(double)RAND_MAX);
#else
    latent_curr[i] = 0.5;
#endif
    error[i] = 0.0;
  }

  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);
#ifdef COMPUTE_RMSE
  rmse = 0;
#endif

  for (int iter = 0; iter < numIter; iter++){
#ifdef COMPUTE_RMSE
    cout << "sum of squared error: " << rmse/2 << " for iter: " << iter << endl;
    rmse = 0;
#endif
    //edgemap to accumulate error for each node
    edgeMap(GA, Frontier, GD_F<vertex>(latent_curr,error,GA.V),GA.m/20);
    //vertexmap to update the latent vectors
    vertexMap(Frontier,GD_Vertex_F(latent_curr,error,step,lambda));

#ifdef DEBUG2
    double latent_sum = 0.0;
    for(long i=0;i<K*n;i++) latent_sum += latent_curr[i];
    cout << "latent sum: " << latent_sum << endl;
#endif

  }


  Frontier.del();

}
