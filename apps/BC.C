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
#include <vector>

typedef double fType;

int singlethreadedUpdateCountBC_F = 0;
int singlethreadedUpdateCountBC_BACK_F = 0;

int singlethreadedAtomicUpdateCountBC_F = 0;
int singlethreadedAtomicUpdateCountBC_BACK_F = 0;

struct BC_F {
  fType* NumPaths;
  bool* Visited;
  unsigned int* ActiveBitVector;

  BC_F(fType* _NumPaths, bool* _Visited, unsigned int* _ActiveBitVector) : 
    NumPaths(_NumPaths), Visited(_Visited), ActiveBitVector(_ActiveBitVector) {}
  __forceinline bool update(uintE s, uintE d){ //Update function for forward phase
    bool result;

    if(!Check_Bit(ActiveBitVector, s)){
      result = false;
    } else {
#ifdef DEBUG
    singlethreadedUpdateCountBC_F++;
#endif
    
    fType oldV = NumPaths[d];
    NumPaths[d] += NumPaths[s];
    result =( oldV == 0.0);
    }
    return result;
  }
  __forceinline bool updateAtomic (uintE s, uintE d) { //atomic Update, basically an add
    volatile fType oldV, newV; 
#ifdef DEBUG
    singlethreadedAtomicUpdateCountBC_F++;
#endif

    do { 
      oldV = NumPaths[d]; newV = oldV + NumPaths[s];
    } while(!CAS(&NumPaths[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

struct BC_Back_F {
  fType* Dependencies;
  bool* Visited;
  unsigned int* ActiveBitVector;

  BC_Back_F(fType* _Dependencies, bool* _Visited, unsigned int* _ActiveBitVector) : 
    Dependencies(_Dependencies), Visited(_Visited), ActiveBitVector(_ActiveBitVector) {}
  __forceinline bool update(uintE s, uintE d){ //Update function for backwards phase

    bool result;
    
    if(!Check_Bit(ActiveBitVector, s)){
      result =  false;
    } else {
#ifdef DEBUG
      singlethreadedUpdateCountBC_BACK_F++;
#endif
      fType oldV = Dependencies[d];
      Dependencies[d] += Dependencies[s];
      result = ( oldV == 0.0);
    }
    return result;
  }
  
  __forceinline bool updateAtomic (uintE s, uintE d) { //atomic Update
#ifdef DEBUG
    singlethreadedAtomicUpdateCountBC_BACK_F++;
#endif
    volatile fType oldV, newV;
    do {
      oldV = Dependencies[d];
      newV = oldV + Dependencies[s];
    } while(!CAS(&Dependencies[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

//vertex map function to mark visited vertexSubset
struct BC_Vertex_F {
  bool* Visited;
  BC_Vertex_F(bool* _Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

//vertex map function (used on backwards phase) to mark visited vertexSubset
//and add to Dependencies score
struct BC_Back_Vertex_F {
  bool* Visited;
  fType* Dependencies, *inverseNumPaths;
  BC_Back_Vertex_F(bool* _Visited, fType* _Dependencies, fType* _inverseNumPaths) : 
    Visited(_Visited), Dependencies(_Dependencies), inverseNumPaths(_inverseNumPaths) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    Dependencies[i] += inverseNumPaths[i];
    return 1; }};

void SetBitVectorWithFrontier(unsigned int* bit_vector, int bit_vector_length, vertexSubset Frontier, int numVertices){
  // reset all bits of the bit vector
  parallel_for(int i = 0; i < bit_vector_length; i++)
    bit_vector[i] = 0;
  
  if (Frontier.isDense) {
    long m = Frontier.numNonzeros();
#ifdef DEBUG
    cout << "dense frontier: " << endl;
    cout << "m: " << m << endl;
#endif
    bool * active = Frontier.d;
    parallel_for(int i = 0; i < numVertices; i+=32){
      int start = i;
      int end = (((i + 32) < numVertices)? (i+32):numVertices);
      for(int j = start; j < end; j++){
	if (active[j])
	  Set_Bit_Only(bit_vector, j);
      }
    }
  }else {
    long m = Frontier.numNonzeros();
#ifdef DEBUG
    cout << "sparse frontier: " << endl;
    cout << "m: " << m << endl;
#endif
    for(int i = 0; i < m; i++){
      Set_Bit_Only(bit_vector, Frontier.s[i]);
    
     }  
    } 
  }

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r",0);
  long n = GA.n, threshold = GA.m/20;
  fType* NumPaths = newA(fType,n);
  {parallel_for(long i=0;i<n;i++) NumPaths[i] = 0.0;}
  NumPaths[start] = 1.0;

#ifdef DEBUG
    singlethreadedUpdateCountBC_F = 0;
    singlethreadedUpdateCountBC_BACK_F = 0;
    singlethreadedAtomicUpdateCountBC_F = 0;
    singlethreadedAtomicUpdateCountBC_BACK_F = 0;
#endif

  bool* Visited = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
  Visited[start] = 1;
  vertexSubset Frontier(n,start);
 
  vector<vertexSubset> Levels;
  Levels.push_back(Frontier);

  long round = 0;

  #ifdef DEBUG2
  cout << "start: " << start << endl;
  cout << "start in degree: " << GA.V[start].getInDegree() << endl;
  cout << "start out degree: " << GA.V[start].getOutDegree() << endl;
  #endif

  int bit_vector_length = n / 32;
  if(n % 32 != 0)
    bit_vector_length++;
  unsigned int* bitVector = (unsigned int*) _mm_malloc(sizeof(unsigned int) * bit_vector_length,64);

  

  while(!Frontier.isEmpty()){ //first phase
#ifdef DEBUG
    cout << "iter: " << round << endl;
    cout << "numActive: " << Frontier.numNonzeros() << endl;
#endif
    round++;
    
    //startTime();
    SetBitVectorWithFrontier(bitVector, bit_vector_length, Frontier, n);
    //nextTime("Set Bit Vector time");

    vertexSubset output = edgeMap(GA, Frontier, BC_F(NumPaths,Visited, bitVector),threshold);
    vertexMap(output, BC_Vertex_F(Visited)); //mark visited
    Levels.push_back(output); //save frontier onto Levels
    Frontier = output;
  }

  fType* Dependencies = newA(fType,n);
  {parallel_for(long i=0;i<n;i++) Dependencies[i] = 0.0;}

  //invert numpaths
  fType* inverseNumPaths = NumPaths;
  {parallel_for(long i=0;i<n;i++) inverseNumPaths[i] = 1/inverseNumPaths[i];}

  Levels[round].del();
  //reuse Visited
  {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
  Frontier = Levels[round-1];
  vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

  //tranpose graph
  GA.transpose();
 
  for(long r=round-2;r>=0;r--) { //backwards phase
    
    //startTime();
    SetBitVectorWithFrontier(bitVector, bit_vector_length, Frontier, n);
    //nextTime("Set bit vector time");


    vertexSubset output = edgeMap(GA, Frontier, 
				  BC_Back_F(Dependencies,Visited, bitVector),threshold);
    output.del(); Frontier.del();
    Frontier = Levels[r]; //gets frontier from Levels array
    //vertex map to mark visited and update Dependencies scores
    vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
  }
  
  Frontier.del();

  //Update dependencies scores
  parallel_for(long i=0;i<n;i++) {
    Dependencies[i]=(Dependencies[i]-inverseNumPaths[i])/inverseNumPaths[i];
  }
  
#ifdef DEBUG
  cout << Dependencies[0] << endl;
  cout << Dependencies[4194303] << endl;
  cout << " single threaded BC_F non atomic update count: " << singlethreadedUpdateCountBC_F << endl;
  cout << " single threaded BC_BACK_F non atomic update count: " << singlethreadedUpdateCountBC_BACK_F << endl;
  cout << " single threaded BC_F atomic update count: " << singlethreadedAtomicUpdateCountBC_F << endl;
  cout << " single threaded BC_BACK_F atomic update count: " << singlethreadedAtomicUpdateCountBC_BACK_F << endl;
#endif
  
  _mm_free(bitVector);
  free(inverseNumPaths);
  free(Visited);
  free(Dependencies);
}
