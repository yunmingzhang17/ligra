// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"
//#define DEBUG
//#define DEBUG1

#ifdef DEBUG
int updateCount = 0;
int atomicUpdateCount = 0;
#endif

struct BFS_F {
  uintE* Parents;
  BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update (uintE s, uintE d) { //Update
    #ifdef DEBUG
    updateCount++;
    #endif
    if(Parents[d] == UINT_E_MAX) { Parents[d] = s; return 1; }
    else return 0;
  }
  inline bool updateAtomic (uintE s, uintE d){ //atomic version of Update
    #ifdef DEBUG
    atomicUpdateCount++;
    #endif
    return (CAS(&Parents[d],UINT_E_MAX,s));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return (Parents[d] == UINT_E_MAX); } 
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r",0);
  long n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  uintE* Parents = newA(uintE,n);
  parallel_for(long i=0;i<n;i++) Parents[i] = UINT_E_MAX;
  Parents[start] = start;
  vertexSubset Frontier(n,start); //creates initial frontier
  int iter = 0;
  #ifdef DEBUG2
  cout << "start: " << start << endl;
  cout << "start in degree: " << GA.V[start].getInDegree() << endl;
  cout << "start out degree: " << GA.V[start].getOutDegree() << endl;
  #endif

  #ifdef DEBUG
  updateCount = 0;
  atomicUpdateCount = 0;
  #endif

  #ifdef BENCHMARK
  startTime();
  #endif

  while(!Frontier.isEmpty()){ //loop until frontier is empty
    #ifdef DEBUG
    cout << "iter: " << iter << endl;
    cout << "numActive: " << Frontier.numNonzeros() << endl;
    #endif
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents),GA.m/20);    
    Frontier.del();
    Frontier = output; //set new frontier
  }

#ifdef BENCHMARK
  nextTimePerIter("time_for_all_iterations,", 1);
#endif

  #ifdef DEBUG
  cout << "non atomic update count (single threaded test): " << updateCount << endl;
  cout << "atomic update count (single threaded test): " << atomicUpdateCount << endl;
  #endif
 
  Frontier.del();
  free(Parents); 
}
