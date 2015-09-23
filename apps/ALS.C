#define WEIGHTED 1
#include "ligra.h"
#include "eigen_wrapper.hpp"

template <class vertex>
struct ALS_Vertex_F {
  vertex* V;
  vector<float> * latent_factors_inmem;

  ALS_Vertex_F(vertex* _V) :
    V(_V){};
  
  // inline bool update(uintE s, uintE d, intE rating){
  //   cout << "s: " << s << " d: " << d << " r: " << rating << endl;
  // }

  // inline bool updateAtomic(uintE s, uintE d, intE rating){
  //   cout << "s: " << s << " d: " << d << " r: " << rating << endl;
  // }
  // inline bool cond (uintE d) { return cond_true(d);}

  inline bool operator() (uintE i){
    cout << "vertex: " << i << endl;
    uintE d = V[i].getInDegree();
    for (uintE j = 0; j < d; j++){
      uintE ngh = V[i].getInNeighbor(j);
      intE rating = V[i].getInWeight(j);
      cout << "edge src: " << i << " dst: " << ngh << " rating: " << rating << endl; 
    }
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
