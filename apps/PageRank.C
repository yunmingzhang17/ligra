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
#include "math.h"
#define HACK
#define NONORM

template <class vertex>
struct PR_F {
  double* p_curr, *p_next;
  vertex* V;
  PR_F(double* _p_curr, double* _p_next, vertex* _V) : 
    p_curr(_p_curr), p_next(_p_next), V(_V) {}
  inline bool update(uintE s, uintE d){ //update function applies PageRank equation
#ifdef DEBUG1
    cout << "update s: " << s << " d: " << d << endl;
#endif
    p_next[d] += p_curr[s]/V[s].getOutDegree();
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  double damping;
  double addedConstant;
  double* p_curr;
  double* p_next;
  PR_Vertex_F(double* _p_curr, double* _p_next, double _damping, intE n) :
    p_curr(_p_curr), p_next(_p_next), 
#ifndef NONORM
    damping(_damping), addedConstant((1-_damping)*(1/(double)n)){}
#else
  damping(_damping), addedConstant((1-_damping)){}
#endif

  inline bool operator () (uintE i) {
 #ifdef DEBUG1
    cout << "PR_VERTEX_F  i: " << i << endl;
#endif
    p_next[i] = damping*p_next[i] + addedConstant;
    return 1;
  }
};

//resets p
struct PR_Vertex_Reset {
  double* p_curr;
  PR_Vertex_Reset(double* _p_curr) :
    p_curr(_p_curr) {}
  inline bool operator () (uintE i) {
    p_curr[i] = 0.0;
    return 1;
  }
};

//#define DEBUG2

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long maxIters = P.getOptionLongValue("-maxiters",5);
  const intE n = GA.n;
  const double damping = 0.85, epsilon = 0.0000001;
  //printf("graph number of vertices: %d number of edges: %d \n", n, GA.m);

  double one_over_n = 1/(double)n;
#ifdef WORKINGSET_EXPAND
  double* p_curr = newA(double,4*n);
  double* p_next = newA(double, 4*n);
#else
  double* p_curr = newA(double,n);
  double* p_next = newA(double,n);
#endif

#ifndef NONORM
  {parallel_for(long i=0;i<n;i++) p_curr[i] = one_over_n;}//use 1 instead of 1/n
#else
  {parallel_for(long i=0;i<n;i++) p_curr[i] = 1.0;}//use 1 instead of 1/n
#endif

 
  {parallel_for(long i=0;i<n;i++) p_next[i] = 0;} //0 if unchanged
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}
  bool* vertexFrontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) vertexFrontier[i] = 1;}
  
  vertexSubset Frontier(n,n,frontier);
  vertexSubset VertexFrontier(n,n,vertexFrontier);
  double rankSum;

  long iter = 0;
  
  #ifdef BENCHMARK
  startTime();
  #endif

  while(iter++ < maxIters){

    #ifdef DEBUG1  
    for(int i = 0; i < n; i++) {
      printf("iter %d curr: %f next: %f frontier: %d \n ", iter, p_curr[i], p_next[i], frontier[i]);    
    }
    #endif

#ifdef PUSH
    vertexSubset output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V),GA.m/20,DENSE_FORWARD);
#else
    vertexSubset output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.V),GA.m/20);
#endif

#ifdef DEBUG1     
    for(int i = 0; i < n; i++) {
      printf("iter %d curr: %f next: %f frontier: %d \n ", iter, p_curr[i], p_next[i], frontier[i]);    
    }
#endif

    vertexMap(Frontier,PR_Vertex_F(p_curr,p_next,damping,n));
    //compute L1-norm between p_curr and p_next
 #ifdef DEBUG2
    rankSum = sequence::plusReduce(p_curr,n);   
    //cout << "rank sum: " << rankSum << endl;
 #endif
    {parallel_for(long i=0;i<n;i++) {
      p_curr[i] = fabs(p_curr[i]-p_next[i]);
      }}
    double L1_norm = sequence::plusReduce(p_curr,n);
    
    if(L1_norm < epsilon) break;
    //reset p_curr
    vertexMap(Frontier,PR_Vertex_Reset(p_curr));
    swap(p_curr,p_next);
#ifdef HACK
    //Frontier.del(); 
    //Frontier = output;
    output.del();
#else
    Frontier.del(); 
    Frontier = output;
#endif

  }

#ifdef BENCHMARK
  nextTimePerIter("time_per_iter,", iter);
#else
  printf("PageRank took %d iterations, max iter: %d, rank sum: %f\n", iter, maxIters, rankSum);
#endif

  Frontier.del(); free(p_curr); free(p_next); 
}
