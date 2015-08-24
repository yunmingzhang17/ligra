#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

struct symmetricVertex {
#ifndef WEIGHTED
  uintE* neighbors;
#else 
  intE* neighbors;
#endif
  uintT degree;
  void del() {free(neighbors); }
#ifndef WEIGHTED
symmetricVertex(uintE* n, uintT d) 
#else 
symmetricVertex(intE* n, uintT d) 
#endif
: neighbors(n), degree(d) {}
#ifndef WEIGHTED
  uintE getInNeighbor(uintT j) { return neighbors[j]; }
  uintE getOutNeighbor(uintT j) { return neighbors[j]; }
  void setInNeighbors(uintE* _i) { neighbors = _i; }
  void setOutNeighbors(uintE* _i) { neighbors = _i; }
  uintE* getInNeighbors() { return neighbors; }
#else
  //weights are stored in the entry after the neighbor ID
  //so size of neighbor list is twice the degree
  intE getInNeighbor(intT j) { return neighbors[2*j]; }
  intE getOutNeighbor(intT j) { return neighbors[2*j]; }
  intE getInWeight(intT j) { return neighbors[2*j+1]; }
  intE getOutWeight(intT j) { return neighbors[2*j+1]; }
  void setInNeighbors(intE* _i) { neighbors = _i; }
  void setOutNeighbors(intE* _i) { neighbors = _i; }
  intE* getInNeighbors() { return neighbors; }
#endif

  uintT getInDegree() { return degree; }
  uintT getOutDegree() { return degree; }
  void setInDegree(uintT _d) { degree = _d; }
  void setOutDegree(uintT _d) { degree = _d; }
  void flipEdges() {}
};

struct asymmetricVertex {
#ifndef WEIGHTED
  uintE* inNeighbors, *outNeighbors;
#else
  intE* inNeighbors, *outNeighbors;
#endif
  uintT outDegree;
  uintT inDegree;
  void del() {free(inNeighbors); free(outNeighbors);}
#ifndef WEIGHTED
asymmetricVertex(uintE* iN, uintE* oN, uintT id, uintT od) 
#else
asymmetricVertex(intE* iN, intE* oN, uintT id, uintT od) 
#endif
: inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {}
  
#ifndef WEIGHTED
  uintE getInNeighbor(uintT j) { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) { return outNeighbors[j]; }
  void setInNeighbors(uintE* _i) { inNeighbors = _i; }
  void setOutNeighbors(uintE* _i) { outNeighbors = _i; }
  uintE* getInNeighbors() { return inNeighbors; }
#else 
  intE getInNeighbor(uintT j) { return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) { return outNeighbors[2*j]; }
  intE getInWeight(uintT j) { return inNeighbors[2*j+1]; }
  intE getOutWeight(uintT j) { return outNeighbors[2*j+1]; }
  void setInNeighbors(intE* _i) { inNeighbors = _i; }
  void setOutNeighbors(intE* _i) { outNeighbors = _i; }
  intE* getInNeighbors() { return inNeighbors; }
#endif
  uintT getInDegree() { return inDegree; }
  uintT getOutDegree() { return outDegree; }
  void setInDegree(uintT _d) { inDegree = _d; }
  void setOutDegree(uintT _d) { outDegree = _d; }
  void flipEdges() { swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }
};

template <class vertex>
struct graph {
  vertex *V;
  long n;
  long m;
#ifndef WEIGHTED
  uintE* allocatedInplace, * inEdges;
#else
  intE* allocatedInplace, * inEdges;
#endif
  uintE* flags;
  bool transposed;
  graph(vertex* VV, long nn, long mm) 
  : V(VV), n(nn), m(mm), allocatedInplace(NULL), flags(NULL), transposed(false) {}
#ifndef WEIGHTED
  graph(vertex* VV, long nn, long mm, uintE* ai, uintE* _inEdges = NULL) 
#else
  graph(vertex* VV, long nn, long mm, intE* ai, intE* _inEdges = NULL) 
#endif
  : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges), flags(NULL), transposed(false) {}
  void del() {
    if (flags != NULL) free(flags);
    if (allocatedInplace == NULL) 
      for (long i=0; i < n; i++) V[i].del();
    else free(allocatedInplace);
    free(V);
    if(inEdges != NULL) free(inEdges);
  }
  void transpose() {
    if(sizeof(vertex) == sizeof(asymmetricVertex)) {
      parallel_for(long i=0;i<n;i++) {
	V[i].flipEdges();
      }
      transposed = !transposed;
    } 
  }
};
