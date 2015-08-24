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
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
using namespace std;

//#define DEBUG1
#define PROFILE

//*****START FRAMEWORK*****

//*****VERTEX OBJECT*****
struct vertexSubset {
  long n, m;
  uintE* s;
  bool* d;
  bool isDense;

  // make a singleton vertex in range of n
vertexSubset(long _n, intE v) 
: n(_n), m(1), d(NULL), isDense(0) {
  s = newA(uintE,1);
  s[0] = v;
}
  
  //empty vertex set
vertexSubset(long _n) : n(_n), m(0), d(NULL), s(NULL), isDense(0) {}
  // make vertexSubset from array of vertex indices
  // n is range, and m is size of array
vertexSubset(long _n, long _m, uintE* indices) 
: n(_n), m(_m), s(indices), d(NULL), isDense(0) {}
  // make vertexSubset from boolean array, where n is range
vertexSubset(long _n, bool* bits) 
: n(_n), d(bits), s(NULL), isDense(1)  {
  m = sequence::sum(bits,_n); }
  // make vertexSubset from boolean array giving number of true values
vertexSubset(long _n, long _m, bool* bits) 
: n(_n), m(_m), s(NULL), d(bits), isDense(1)  {}

  // delete the contents
  void del(){
    if (d != NULL) free(d);
    if (s != NULL) free(s);
  }
  long numRows() { return n; }
  long numNonzeros() { return m; }
  bool isEmpty() { return m==0; }

  // converts to dense but keeps sparse representation if there
  void toDense() {
    if (d == NULL) {
      d = newA(bool,n);
      {parallel_for(long i=0;i<n;i++) d[i] = 0;}
      {parallel_for(long i=0;i<m;i++) d[s[i]] = 1;}
    }
    isDense = true;
  }

  // converts to sparse but keeps dense representation if there
  void toSparse() {
    if (s == NULL) {
      _seq<uintE> R = sequence::packIndex<uintE>(d,n);
      if (m != R.n) {
	cout << "bad stored value of m" << endl;
	cout << "m: " << m <<  endl;
	abort();
      }
      s = R.A;
    }
    isDense = false;
  }
  // check for equality
  bool eq (vertexSubset& b) {
    toDense();
    b.toDense();
    bool* c = newA(bool,n);
    {parallel_for (long i=0; i<b.n; i++) 
	c[i] = (d[i] != b.d[i]);}
    bool equal = (sequence::sum(c,n) == 0);
    free(c);
    return equal;
  }

  void print() {
    if (isDense) {
      cout << "D:";
      for (long i=0;i<n;i++) if (d[i]) cout << i << " ";
      cout << endl;
    } else {
      cout << "S:";
      for (long i=0; i<m; i++)	cout << s[i] << " ";
      cout << endl;
    }
  }
};

struct nonMaxF{bool operator() (uintE &a) {return (a != UINT_E_MAX);}};

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options { DENSE, DENSE_FORWARD};

//remove duplicate integers in [0,...,n-1]
void remDuplicates(uintE* indices, uintE* flags, long m, long n) {
  //make flags for first time
  if(flags == NULL) {flags = newA(uintE,n); 
    {parallel_for(long i=0;i<n;i++) flags[i]=UINT_E_MAX;}}
  {parallel_for(uintE i=0;i<m;i++)
      if(indices[i] != UINT_E_MAX && flags[indices[i]] == UINT_E_MAX) 
	CAS(&flags[indices[i]],(uintE)UINT_E_MAX,i);
  }
  //reset flags
  {parallel_for(long i=0;i<m;i++){
      if(indices[i] != UINT_E_MAX){
	if(flags[indices[i]] == i){ //win
	  flags[indices[i]] = UINT_E_MAX; //reset
	}
	else indices[i] = UINT_E_MAX; //lost
      }
    }
  }
}

//*****EDGE FUNCTIONS*****
template <class F, class vertex>
  bool* recursiveEdgeMapDense(graph<vertex> GA, bool* vertexSubset, F f, bool parallel = 0) {
  cout << "recursive " << endl;
#ifdef DEBUG1
  cout << "inside recrusiveEdgeMapDense" << endl;
#endif
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  cilk_spawn recursiveEdgeMapDenseHelper(GA, 0, numVertices-1, vertexSubset, f, next);


  //process the last node
  for(int j = 0; j < GA.V[numVertices-1].getInDegree(); j++){
    next[numVertices-1] = 0;
    uintE ngh = GA.V[numVertices-1].getInNeighbor(j);
#ifdef DEBUG1
    cout << "node: " << (numVertices-1) << " ngh: " << ngh << endl;
#endif
#ifndef WEIGHTED
    if (vertexSubset[ngh] && f.update(ngh,(numVertices-1)))
#else
    if (vertexSubset[ngh] && f.update(ngh,(numVertices-1),GA.V[(numVertices-1)].getInWeight(j)))
#endif
      next[(numVertices-1)] = 1;
    if(!f.cond((numVertices-1))) break;
  }

  cilk_sync;
  return next;
}

template <class F, class vertex>
  void recursiveEdgeMapDenseHelper(graph<vertex> GA, int start, int end, bool* vertexSubset, F f, bool *next, int edgeGrainSize = 4096){
  //int edgeEnd = reinterpret_cast<int>(GA.V[end].getInNeighbors());
  //int edgeStart = reinterpret_cast<int>(GA.V[start].getInNeighbors());
  //cout << GA.V[end].getInNeighbors() << endl;
  //cout << GA.V[start].getInNeighbors() << endl;
#ifdef DEBUG1
  cout << "start: " << start << " end: " << end << endl;
#endif

  long diff = reinterpret_cast<long>((GA.V[end].getInNeighbors() - GA.V[start].getInNeighbors()));
#ifdef DEBUG1
  cout << "diff: " << diff << endl;
#endif
  if ((start == end -1) || (diff < edgeGrainSize)){
    for (int i = start; i < end; i++){
      next[i] = 0;
      for(int j = 0; j < GA.V[i].getInDegree(); j++){
	uintE ngh = GA.V[i].getInNeighbor(j);
#ifdef DEBUG1
	cout << "node: " << i << " ngh: " << ngh << endl;
#endif

#ifndef WEIGHTED
	if (vertexSubset[ngh] && f.update(ngh,i))
#else
	if (vertexSubset[ngh] && f.update(ngh,i,GA.V[i].getInWeight(j)))
#endif
	  next[i] = 1;
	if(!f.cond(i)) break;
      }
    }
  } else {
	cilk_spawn recursiveEdgeMapDenseHelper(GA, start, start + ((end-start)>>1), vertexSubset, f, next, edgeGrainSize);
	recursiveEdgeMapDenseHelper(GA, start + ((end-start)>>1), end, vertexSubset, f, next, edgeGrainSize);
  }
}



template <class F, class vertex>
  bool* edgeMapDense(graph<vertex> GA, bool* vertexSubset, F f, bool parallel = 0) {
#ifdef DEBUG1
  cout << "inside EdgeMapDense" << endl;
#endif
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for (long i=0; i<numVertices; i++){
    next[i] = 0;
    if (f.cond(i)) { 
      uintE d = G[i].getInDegree();
      if(!parallel || d < 1000) {
	for(uintE j=0; j<d; j++){
	  uintE ngh = G[i].getInNeighbor(j);
#ifdef DEBUG1
	  cout << "node: " << i <<  " ngh: " << ngh << endl;
#endif

#ifndef WEIGHTED
	  if (vertexSubset[ngh] && f.update(ngh,i))
#else
	  if (vertexSubset[ngh] && f.update(ngh,i,G[i].getInWeight(j)))
#endif
	    next[i] = 1;
	  if(!f.cond(i)) break;
	}
      } else {
	{parallel_for(uintE j=0; j<d; j++){
	  uintE ngh = G[i].getInNeighbor(j);
#ifndef WEIGHTED
	  if (vertexSubset[ngh] && f.update(ngh,i))
#else
	  if (vertexSubset[ngh] && f.update(ngh,i,G[i].getInWeight(j)))
#endif
	    next[i] = 1;
	  }}
      }
    }
    }}

#ifdef DEBUG1
  cout << "next " << endl;
  for (int i = 0; i < numVertices; i++) {
    cout << next[i] << " " ;
  }
  cout << endl;
#endif
  return next;
}

template <class F, class vertex>
bool* edgeMapDenseForward(graph<vertex> GA, bool* vertexSubset, F f) {
#ifdef DEBUG1
  cout << "inside EdgeMapDenseForward" << endl;
#endif
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for(long i=0;i<numVertices;i++) next[i] = 0;}
  {parallel_for (long i=0; i<numVertices; i++){
    if (vertexSubset[i]) {
      uintE d = G[i].getOutDegree();
      if(d < 1000) {
	for(uintE j=0; j<d; j++){
	  uintE ngh = G[i].getOutNeighbor(j);
#ifndef WEIGHTED
	  if (f.cond(ngh) && f.updateAtomic(i,ngh))
#else 
	  if (f.cond(ngh) && f.updateAtomic(i,ngh,G[i].getOutWeight(j))) 
#endif
	    next[ngh] = 1;
	}
      }
      else {
	{parallel_for(uintE j=0; j<d; j++){
	  uintE ngh = G[i].getOutNeighbor(j);
#ifndef WEIGHTED
	  if (f.cond(ngh) && f.updateAtomic(i,ngh)) 
#else
	    if (f.cond(ngh) && f.updateAtomic(i,ngh,G[i].getOutWeight(j)))
#endif
	  next[ngh] = 1;
	  }}
      }
    }
    }}
  return next;
}

template <class F, class vertex>
pair<long,uintE*> edgeMapSparse(vertex* frontierVertices, uintE* indices, 
				uintT* degrees, uintT m, F f, 
				long remDups=0, uintE* flags=NULL) {
#ifdef DEBUG1
  cout << "inside EdgeMapSparse" << endl;
#endif
  uintT* offsets = degrees;
  long outEdgeCount = sequence::plusScan(offsets, degrees, m);
  uintE* outEdges = newA(uintE,outEdgeCount);
  {parallel_for (long i = 0; i < m; i++) {
      uintT v = indices[i], o = offsets[i];
    vertex vert = frontierVertices[i]; 
    uintE d = vert.getOutDegree();
    if(d < 1000) {
      for (uintE j=0; j < d; j++) {
	uintE ngh = vert.getOutNeighbor(j);
#ifndef WEIGHTED
	if(f.cond(ngh) && f.updateAtomic(v,ngh)) 
#else
	if(f.cond(ngh) && f.updateAtomic(v,ngh,vert.getOutWeight(j)))
#endif
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = UINT_E_MAX;
      } 
    } else {
      {parallel_for (uintE j=0; j < d; j++) {
	uintE ngh = vert.getOutNeighbor(j);
#ifndef WEIGHTED
	if(f.cond(ngh) && f.updateAtomic(v,ngh)) 
#else
	if(f.cond(ngh) && f.updateAtomic(v,ngh,vert.getOutWeight(j)))
#endif
	  outEdges[o+j] = ngh;
	else outEdges[o+j] = UINT_E_MAX;
	}} 
    }
    }}
  uintE* nextIndices = newA(uintE, outEdgeCount);
  if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);
  // Filter out the empty slots (marked with -1)
  long nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonMaxF());
  free(outEdges);
  return pair<long,uintE*>(nextM, nextIndices);
}

static int edgesTraversed = 0;

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
vertexSubset edgeMap(graph<vertex> GA, vertexSubset &V, F f, intT threshold = -1, 
		 char option=DENSE, bool remDups=false) {
  long numVertices = GA.n, numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  long m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }

  #ifdef DEBUG1
  cout << "Inside edgeMap" << endl;
  #endif

  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(vertex,m);
  {parallel_for (long i=0; i < m; i++){
    vertex v = G[V.s[i]];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  uintT outDegrees = sequence::plusReduce(degrees, m);
  edgesTraversed += outDegrees;
  if (outDegrees == 0) return vertexSubset(numVertices);
  if (m + outDegrees > threshold) { 
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R;
    if (option == DENSE_FORWARD) { 
	R = edgeMapDenseForward(GA,V.d,f); 
	} else { 
#ifndef RECURSIVE
	R = edgeMapDense(GA, V.d, f, option);
#else
	R = recursiveEdgeMapDense(GA, V.d, f, option);
#endif
      }
    vertexSubset v1 = vertexSubset(numVertices, R);
    //cout << "size (D) = " << v1.m << endl;
    return  v1;
  } else { 
    pair<long,uintE*> R = 
      remDups ? 
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    //cout << "size (S) = " << R.first << endl;
    free(degrees);
    free(frontierVertices);
    return vertexSubset(numVertices, R.first, R.second);
  }
}

//*****VERTEX FUNCTIONS*****

//Note: this is the optimized version of vertexMap which does not
//perform a filter
template <class F>
void vertexMap(vertexSubset V, F add) {
  long n = V.numRows(), m = V.numNonzeros();
  if(V.isDense) {
    {parallel_for(long i=0;i<n;i++)
	if(V.d[i]) add(i);}
  } else {
    {parallel_for(long i=0;i<m;i++)
	add(V.s[i]);}
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine);

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool binary = P.getOptionValue("-b");
  long rounds = P.getOptionLongValue("-rounds",3);
  if(symmetric) {
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,symmetric,binary); //symmetric graph
    Compute(G,P);
    for(int r=0;r<rounds;r++) {
      startTime();
      Compute(G,P);
      nextTime("Running time");
    }
    G.del();
  } else {
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,symmetric,binary); //asymmetric graph
    Compute(G,P);
    if(G.transposed) G.transpose();
#ifdef PROFILE
    ofstream outputFile;
    outputFile.open("./.ready");
    outputFile.close();
    cout << "continue?" << "\n";
    int c = getchar();
#endif
    for(int r=0;r<rounds;r++) {
      startTime();
      Compute(G,P);
      nextTime("Running time");
      if(G.transposed) G.transpose();
    }
    G.del();
  }
}
#endif
