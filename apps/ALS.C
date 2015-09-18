#define WEIGHTED 1
#include "ligra.h"

struct ALS_Vertex_F {
  ALS_Vertex_F() {};
  
  // inline bool update(uintE s, uintE d, intE rating){
  //   cout << "s: " << s << " d: " << d << " r: " << rating << endl;
  // }

  // inline bool updateAtomic(uintE s, uintE d, intE rating){
  //   cout << "s: " << s << " d: " << d << " r: " << rating << endl;
  // }
  // inline bool cond (uintE d) { return cond_true(d);}

  inline bool operator() (uintE i){
    cout << "vertex: " << i << endl;
  }
};
  
template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long n = GA.n;
  cout << "num vertices: " << n << endl;
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  vertexSubset Frontier(n,n,frontier);
  vertexMap(Frontier, ALS_Vertex_F());
  Frontier.del();
}
