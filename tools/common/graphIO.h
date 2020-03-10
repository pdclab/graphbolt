// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
//
// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#ifndef _BENCH_GRAPH_IO
#define _BENCH_GRAPH_IO

#include "../../core/common/blockRadixSort.h"
#include "../../core/common/parallel.h"
#include "../../core/common/quickSort.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdint.h>
using namespace std;

typedef pair<long, long> intPair;
typedef pair<long, pair<long, long>> intTriple;

// **************************************************************
//    EDGE ARRAY REPRESENTATION
// **************************************************************

template <class intT> struct edge {
  intT u;
  intT v;
  edge(intT f, intT s) : u(f), v(s) {}
  edge() : u(0), v(0) {}
};

template <class intT> struct edgeArray {
  edge<intT> *E;
  intT numRows;
  intT numCols;
  intT nonZeros;
  void del() { free(E); }
  edgeArray(edge<intT> *EE, intT r, intT c, intT nz)
      : E(EE), numRows(r), numCols(c), nonZeros(nz) {}
  edgeArray() {}
};

template <class intT> struct wghEdge {
  intT u;
  intT v;
  intT w;
  wghEdge(intT f, intT s, intT t) : u(f), v(s), w(t) {}
};

template <class intT> struct wghEdgeArray {
  wghEdge<intT> *E;
  intT numRows;
  intT numCols;
  intT nonZeros;
  void del() { free(E); }
  wghEdgeArray(wghEdge<intT> *EE, intT r, intT c, intT nz)
      : E(EE), numRows(r), numCols(c), nonZeros(nz) {}
  wghEdgeArray() {}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intT> struct vertex {
  intT *Neighbors;
  intT degree;
  void del() { free(Neighbors); }
  vertex(intT *N, intT d) : Neighbors(N), degree(d) {}
};

template <class intT> struct wghVertex {
  intT *Neighbors;
  intT degree;
  intT *nghWeights;
  void del() {
    free(Neighbors);
    free(nghWeights);
  }
  wghVertex(intT *N, intT *W, intT d)
      : Neighbors(N), nghWeights(W), degree(d) {}
};

template <class intT> struct graph {
  vertex<intT> *V;
  intT n;
  intT m;
  intT *allocatedInplace;
  graph(vertex<intT> *VV, intT nn, uintT mm)
      : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
  graph(vertex<intT> *VV, intT nn, uintT mm, intT *ai)
      : V(VV), n(nn), m(mm), allocatedInplace(ai) {}
  intT *vertices() { return allocatedInplace + 2; }
  intT *edges() { return allocatedInplace + 2 + n; }
  graph copy() {
    vertex<intT> *VN = newA(vertex<intT>, n);
    intT *_allocatedInplace = newA(intT, n + m + 2);
    _allocatedInplace[0] = n;
    _allocatedInplace[1] = m;
    intT *Edges = _allocatedInplace + n + 2;
    intT k = 0;
    for (intT i = 0; i < n; i++) {
      _allocatedInplace[i + 2] = allocatedInplace[i + 2];
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      for (intT j = 0; j < V[i].degree; j++)
        Edges[k++] = V[i].Neighbors[j];
    }
    return graph(VN, n, m, _allocatedInplace);
  }
  void del() {
    if (allocatedInplace == NULL)
      for (intT i = 0; i < n; i++)
        V[i].del();
    else
      free(allocatedInplace);
    free(V);
  }
};

template <class intT> struct wghGraph {
  wghVertex<intT> *V;
  intT n;
  uintT m;
  intT *allocatedInplace;
  intT *weights;
  wghGraph(wghVertex<intT> *VV, intT nn, uintT mm)
      : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
  wghGraph(wghVertex<intT> *VV, intT nn, uintT mm, intT *ai, intT *_weights)
      : V(VV), n(nn), m(mm), allocatedInplace(ai), weights(_weights) {}
  wghGraph copy() {
    wghVertex<intT> *VN = newA(wghVertex<intT>, n);
    intT *Edges = newA(intT, m);
    intT *Weights = newA(intT, m);
    intT k = 0;
    for (intT i = 0; i < n; i++) {
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      VN[i].nghWeights = Weights + k;
      for (intT j = 0; j < V[i].degree; j++) {
        Edges[k] = V[i].Neighbors[j];
        Weights[k++] = V[i].nghWeights[j];
      }
    }
    return wghGraph(VN, n, m, Edges, Weights);
  }
  void del() {
    if (allocatedInplace == NULL)
      for (intT i = 0; i < n; i++)
        V[i].del();
    else {
      free(allocatedInplace);
    }
    free(V);
  }
};

// **************************************************************
//    GRAPH UTILITIES
// **************************************************************
struct edgeCmp {
  bool operator()(edge<uintT> e1, edge<uintT> e2) {
    return ((e1.u < e2.u) ? 1 : ((e1.u > e2.u) ? 0 : (e1.v < e2.v)));
  }
};

struct edgeCopy {
  void operator()(edge<uintT> e1, edge<uintT> e2) {
    e2.u = e1.u;
    e2.v = e1.v;
  }
};

struct edgeCmpWgh {
  bool operator()(wghEdge<uintT> e1, wghEdge<uintT> e2) {
    return ((e1.u < e2.u) ? 1 : ((e1.u > e2.u) ? 0 : (e1.v < e2.v)));
  }
};

template <class intT> void printEdges(edge<intT> *E, intT m) {
  cout << "Edges:\n";
  for (intT i = 0; i < m; i++) {
    cout << E[i].u << "," << E[i].v << endl;
  }
}

template <class intT> edgeArray<intT> remDuplicates(edgeArray<intT> A) {
  intT m = A.nonZeros;
  edge<intT> *E = newA(edge<intT>, m);
  {
    parallel_for(intT i = 0; i < m; i++) {
      E[i].u = A.E[i].u;
      E[i].v = A.E[i].v;
    }
  }
  quickSort(E, m, edgeCmp());
  intT *flags = newA(intT, m);
  flags[0] = 1;
  {
    parallel_for(intT i = 1; i < m; i++) {
      if ((E[i].u != E[i - 1].u) || (E[i].v != E[i - 1].v))
        flags[i] = 1;
      else
        flags[i] = 0;
    }
  }

  intT mm = sequence::plusScan(flags, flags, m);
  edge<intT> *F = newA(edge<intT>, mm);
  F[mm - 1] = E[m - 1];
  {
    parallel_for(intT i = 0; i < m - 1; i++) {
      if (flags[i] != flags[i + 1])
        F[flags[i]] = E[i];
    }
  }
  free(flags);
  free(E);
  return edgeArray<intT>(F, A.numRows, A.numCols, mm);
}

template <class intT> wghEdgeArray<intT> remDuplicates(wghEdgeArray<intT> A) {
  intT m = A.nonZeros;
  wghEdge<intT> *E = newA(wghEdge<intT>, m);
  {
    parallel_for(intT i = 0; i < m; i++) {
      E[i].u = A.E[i].u;
      E[i].v = A.E[i].v;
      E[i].w = A.E[i].w;
    }
  }
  quickSort(E, m, edgeCmpWgh());

  intT *flags = newA(intT, m);
  flags[0] = 1;
  {
    parallel_for(intT i = 1; i < m; i++) {
      if ((E[i].u != E[i - 1].u) || (E[i].v != E[i - 1].v))
        flags[i] = 1;
      else
        flags[i] = 0;
    }
  }

  intT mm = sequence::plusScan(flags, flags, m);
  wghEdge<intT> *F = newA(wghEdge<intT>, mm);
  F[mm - 1] = E[m - 1];
  {
    parallel_for(intT i = 0; i < m - 1; i++) {
      if (flags[i] != flags[i + 1])
        F[flags[i]] = E[i];
    }
  }
  free(flags);
  return wghEdgeArray<intT>(F, A.numRows, A.numCols, mm);
}

template <class intT> struct nEQF {
  bool operator()(edge<intT> e) { return (e.u != e.v); }
};

template <class intT> struct nEQFWgh {
  bool operator()(wghEdge<intT> e) { return (e.u != e.v); }
};

template <class intT> edgeArray<intT> makeSymmetric(edgeArray<intT> A) {
  intT m = A.nonZeros;
  edge<intT> *E = A.E;
  edge<intT> *F = newA(edge<intT>, 2 * m);
  intT mm = sequence::filter(E, F, m, nEQF<intT>());

  parallel_for(intT i = 0; i < mm; i++) {
    F[i + mm].u = F[i].v;
    F[i + mm].v = F[i].u;
  }

  edgeArray<intT> R =
      remDuplicates(edgeArray<intT>(F, A.numRows, A.numCols, 2 * mm));
  free(F);

  return R;
}

template <class intT> wghEdgeArray<intT> makeSymmetric(wghEdgeArray<intT> A) {
  intT m = A.nonZeros;
  wghEdge<intT> *E = A.E;
  wghEdge<intT> *F = newA(wghEdge<intT>, 2 * m);
  intT mm = sequence::filter(E, F, m, nEQFWgh<intT>());

  parallel_for(intT i = 0; i < mm; i++) {
    F[i + mm].u = F[i].v;
    F[i + mm].v = F[i].u;
    F[i + mm].w = F[i].w;
  }

  wghEdgeArray<intT> R =
      remDuplicates(wghEdgeArray<intT>(F, A.numRows, A.numCols, 2 * mm));
  free(F);

  return R;
}

template <class intT> struct getuF {
  intT operator()(edge<intT> e) { return e.u; }
};

template <class E> struct pairFirstCmp {
  intT operator()(E e) { return e.u; }

  bool operator()(E a, E b) { return a.u < b.u; }
};
// template <class intT>
// struct pairFirstCmp
// {
//   intT operator()(edge<intT> e) { return e.u; }

//   bool operator()(edge<intT> a, edge<intT> b)
//   {
//     return a.u < b.u;
//   }
// };

template <class E> struct pairFirstCmp1 {
  bool operator()(pair<E, E> a, pair<E, E> b) { return a.first < b.first; }
};

template <class intT> struct getuFWgh {
  intT operator()(wghEdge<intT> e) { return e.u; }
};

template <class intT>
graph<intT> graphFromEdges(edgeArray<intT> EA, bool makeSym) {
  edgeArray<intT> A;
  if (makeSym) {
    A = makeSymmetric<intT>(EA);
  } else { // should have copy constructor
    edge<intT> *E = newA(edge<intT>, EA.nonZeros);
    parallel_for(intT i = 0; i < EA.nonZeros; i++) E[i] = EA.E[i];
    A = edgeArray<intT>(E, EA.numRows, EA.numCols, EA.nonZeros);
  }
  intT m = A.nonZeros;
  intT n = max<intT>(A.numCols, A.numRows);
  intT *offsets = newA(intT, n * 2);
  // quickSort(A.E, m, pairFirstCmp<intT>());
  // quickSort(A.E, m, pairFirstCmp<edge<int>>());
  quickSort(A.E, m,
            [](edge<intT> val1, edge<intT> val2) { return val1.u < val2.u; });
  parallel_for(intT i = 0; i < n; i++) { offsets[i] = m; }
  parallel_for(intT i = 0; i < m - 1; i++) {
    long currV = A.E[i].u;
    long nextV = A.E[i + 1].u;
    if (currV != nextV) {
      offsets[nextV] = i + 1;
    }
  }
  offsets[A.E[0].u] = 0;
  sequence::scanIBack(offsets, offsets, n, minF<intT>(), (intT)m);

  // #endif
  intT *X = newA(intT, m);
  vertex<intT> *v = newA(vertex<intT>, n);
  parallel_for(intT i = 0; i < n; i++) {
    intT o = offsets[i];
    intT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = X + o;
    for (intT j = 0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o + j].v;
    }
    if (l > 1) {
      quickSort(v[i].Neighbors, l,
                [](intT val1, intT val2) { return val1 < val2; });
    }
  }
  A.del();
  EA.del();
  // free(offsets);
  return graph<intT>(v, n, m, X);
}

template <class intT>
wghGraph<intT> wghGraphFromWghEdges(wghEdgeArray<intT> EA, bool makeSym) {
  wghEdgeArray<intT> A;
  if (makeSym)
    A = makeSymmetric<intT>(EA);
  else { // should have copy constructor
    wghEdge<intT> *E = newA(wghEdge<intT>, EA.nonZeros);
    parallel_for(intT i = 0; i < EA.nonZeros; i++) E[i] = EA.E[i];
    A = wghEdgeArray<intT>(E, EA.numRows, EA.numCols, EA.nonZeros);
  }
  intT m = A.nonZeros;
  intT n = max<intT>(A.numCols, A.numRows);
  intT *offsets = newA(intT, n * 2);
  intSort::iSort(A.E, offsets, m, n, getuFWgh<intT>());
  intT *X = newA(intT, 2 * m);
  intT *Weights = X + m;
  wghVertex<intT> *v = newA(wghVertex<intT>, n);
  parallel_for(intT i = 0; i < n; i++) {
    intT o = offsets[i];
    intT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = X + o;
    v[i].nghWeights = Weights + o;
    for (intT j = 0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o + j].v;
      v[i].nghWeights[j] = A.E[o + j].w;
    }
  }
  A.del();
  free(offsets);
  return wghGraph<intT>(v, n, m, X, Weights);
}

// **************************************************************
//    BASIC I/O
// **************************************************************
namespace benchIO {
using namespace std;

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n;         // total number of characters
  char *Chars;    // array storing all strings
  long m;         // number of substrings
  char **Strings; // pointers to strings (all should be null terminated)
  words() {}
  words(char *C, long nn, char **S, long mm)
      : Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {
    free(Chars);
    free(Strings);
  }
};

inline bool isSpace(char c) {
  switch (c) {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ':
    return true;
  default:
    return false;
  }
}

struct toLong {
  long operator()(bool v) { return (long)v; }
};

// parallel code for converting a string to words
words stringToWords(char *Str, unsigned long n) {
  parallel_for(unsigned long i = 0; i < n; i++) if (isSpace(Str[i])) Str[i] = 0;

  // mark start of words
  bool *FL = newA(bool, n);
  FL[0] = Str[0];
  parallel_for(unsigned long i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1];

  // offset for each start of word
  _seq<unsigned long> Off = sequence::packIndex<unsigned long>(FL, n);
  unsigned long m = Off.n;
  unsigned long *offsets = Off.A;
  free(FL);

  // pointer to each start of word
  char **SA = newA(char *, m);
  parallel_for(unsigned long j = 0; j < m; j++) SA[j] = Str + offsets[j];

  free(offsets);
  return words(Str, n, SA, m);
}

int writeStringToFile(char *S, unsigned long n, char *fileName) {
  ofstream file(fileName, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return 1;
  }
  file.write(S, n);
  file.close();
  return 0;
}

inline int xToStringLen(long a) { return 21; }
inline void xToString(char *s, long a) { sprintf(s, "%ld", a); }

inline int xToStringLen(unsigned long a) { return 21; }
inline void xToString(char *s, unsigned long a) { sprintf(s, "%lu", a); }

inline uint xToStringLen(uint a) { return 12; }
inline void xToString(char *s, uint a) { sprintf(s, "%u", a); }

inline int xToStringLen(int a) { return 12; }
inline void xToString(char *s, int a) { sprintf(s, "%d", a); }

inline int xToStringLen(double a) { return 18; }
inline void xToString(char *s, double a) { sprintf(s, "%.11le", a); }

inline int xToStringLen(char *a) { return strlen(a) + 1; }
inline void xToString(char *s, char *a) { sprintf(s, "%s", a); }

template <class A, class B> inline int xToStringLen(pair<A, B> a) {
  return xToStringLen(a.first) + xToStringLen(a.second) + 1;
}
template <class A, class B> inline void xToString(char *s, pair<A, B> a) {
  int l = xToStringLen(a.first);
  xToString(s, a.first);
  s[l] = ' ';
  xToString(s + l + 1, a.second);
}

struct notZero {
  bool operator()(char A) { return A > 0; }
};

template <class T> _seq<char> arrayToString(T *A, long n) {
  long *L = newA(long, n);
  { parallel_for(long i = 0; i < n; i++) L[i] = xToStringLen(A[i]) + 1; }
  long m = sequence::scan(L, L, n, addF<long>(), (long)0);
  char *B = newA(char, m);
  parallel_for(long j = 0; j < m; j++) B[j] = 0;
  parallel_for(long i = 0; i < n - 1; i++) {
    xToString(B + L[i], A[i]);
    B[L[i + 1] - 1] = '\n';
  }
  xToString(B + L[n - 1], A[n - 1]);
  B[m - 1] = '\n';
  free(L);
  char *C = newA(char, m + 1);
  long mm = sequence::filter(B, C, m, notZero());
  C[mm] = 0;
  free(B);
  return _seq<char>(C, mm);
}

template <class T> void writeArrayToStream(ofstream &os, T *A, long n) {
  long BSIZE = 1000000;
  long offset = 0;
  while (offset < n) {
    // Generates a string for a sequence of size at most BSIZE
    // and then wrties it to the output stream
    _seq<char> S = arrayToString(A + offset, min(BSIZE, n - offset));
    os.write(S.A, S.n);
    S.del();
    offset += BSIZE;
  }
}

template <class T>
int writeArrayToFile(string header, T *A, long n, char *fileName) {
  ofstream file(fileName, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return 1;
  }
  file << header << endl;
  writeArrayToStream(file, A, n);
  file.close();
  return 0;
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file(fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg(0, ios::beg);
  long n = end - file.tellg();
  char *bytes = newA(char, n + 1);
  file.read(bytes, n);
  file.close();
  return _seq<char>(bytes, n);
}
}; // namespace benchIO

// **************************************************************
//    GRAPH I/O
// **************************************************************
using namespace benchIO;

template <class intT> int xToStringLen(edge<intT> a) {
  return xToStringLen(a.u) + xToStringLen(a.v) + 1;
}

template <class intT> void xToString(char *s, edge<intT> a) {
  int l = xToStringLen(a.u);
  xToString(s, a.u);
  s[l] = ' ';
  xToString(s + l + 1, a.v);
}

namespace benchIO {
using namespace std;

string AdjGraphHeader = "AdjacencyGraph";
string WghAdjGraphHeader = "WeightedAdjacencyGraph";

template <class intT> int writeGraphToFile(graph<intT> G, char *fname) {
  long m = G.m;
  long n = G.n;
  long totalLen = 2 + n + m;
  intT *Out = newA(uintT, totalLen);
  Out[0] = n;
  Out[1] = m;
  parallel_for(long i = 0; i < n; i++) { Out[i + 2] = G.V[i].degree; }
  long total = sequence::scan(Out + 2, Out + 2, n, addF<intT>(), (intT)0);
  for (long i = 0; i < n; i++) {
    intT *O = Out + (2 + n + Out[i + 2]);
    vertex<intT> v = G.V[i];
    for (long j = 0; j < v.degree; j++)
      O[j] = v.Neighbors[j];
  }
  int r = writeArrayToFile(AdjGraphHeader, Out, totalLen, fname);
  free(Out);
  // G.del();
  return r;
}

template <class T> void readBinaryFile(string fileName) {
  std::ifstream iFile1(fileName.c_str(), ios::in | ios::binary);
  iFile1.seekg(0, ios::end);
  T size = iFile1.tellg();
  iFile1.seekg(0);
  cout << "Reading binary file\n";
  cout << "Size : " << size << "\n";
  char *s = (char *)malloc(size);
  iFile1.read(s, size);
  iFile1.close();
  // cout << "s : " << *((T *)s) << "\n";
  T *numbers = (T *)s;
  T len = size / sizeof(T);
  cout << "Numbers : " << len << "\n";
}

template <class intT, class T>
int writeGraphToBinary(graph<intT> G, char *fname) {
  T m = G.m;
  T n = G.n;
  T totalLen = 2 + n + m;
  T *Out = newA(T, totalLen);
  T *outOffsets = &Out[2];
  Out[0] = n;
  Out[1] = m;
  parallel_for(T i = 0; i < n; i++) { Out[i + 2] = (T)G.V[i].degree; }
  T total = sequence::scan(Out + 2, Out + 2, n, addF<T>(), (T)0);
  for (T i = 0; i < n; i++) {
    T *O = Out + (2 + n + Out[i + 2]);
    vertex<intT> v = G.V[i];
    for (T j = 0; j < v.degree; j++)
      O[j] = (T)v.Neighbors[j];
  }

  string fileName(fname);
  cout << "fileName : " << fileName << "\n";

  ofstream file_stream(fileName + ".csr", ios::out | ios::binary);
  if (!file_stream.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return 1;
  }
  int intSize = sizeof(int64_t);

  file_stream.write((char *)Out, sizeof(T) * totalLen);
  file_stream.close();
  readBinaryFile<T>(fileName + ".csr");

  // -----------------------------------

  // Conver csr to csc

  T *inOffsets = newA(T, n);
  parallel_for(T i = 0; i < n; i++) { inOffsets[i] = INT_T_MAX; }
  typedef pair<T, T> longPair;
  longPair *temp = newA(longPair, m);

  parallel_for(T i = 0; i < n; i++) {
    T o = outOffsets[i];
    vertex<intT> v = G.V[i];
    for (T j = 0; j < v.degree; j++) {
      temp[o + j] = make_pair((T)v.Neighbors[j], i);
    }
  }
  quickSort(temp, m, pairFirstCmp1<T>());
  inOffsets[temp[0].first] = 0;
  T *inEdges = newA(T, m);
  inEdges[0] = temp[0].second;
  parallel_for(T i = 1; i < m; i++) {
    inEdges[i] = temp[i].second;
    if (temp[i].first != temp[i - 1].first) {
      inOffsets[temp[i].first] = i;
    }
  }
  sequence::scanIBack(inOffsets, inOffsets, n, minF<T>(), (T)m);

  parallel_for(T i = 0; i < n; i++) { Out[i + 2] = (T)inOffsets[i]; }

  parallel_for(T i = 0; i < m; i++) { Out[i + n + 2] = (T)inEdges[i]; }

  ofstream file_stream2(fileName + ".csc", ios::out | ios::binary);
  if (!file_stream2.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return 1;
  }

  file_stream2.write((char *)Out, sizeof(T) * totalLen);
  file_stream2.close();
  readBinaryFile<T>(fileName + ".csc");

  // G.del();
  // int r = writeArrayToFile(AdjGraphHeader, Out, totalLen, fname);
  free(Out);
  return 0;
}

template <class intT> int writeWghGraphToFile(wghGraph<intT> G, char *fname) {
  long m = G.m;
  long n = G.n;
  long totalLen = 2 + n + m * 2;
  intT *Out = newA(intT, totalLen);
  Out[0] = n;
  Out[1] = m;
  parallel_for(long i = 0; i < n; i++) { Out[i + 2] = G.V[i].degree; }
  long total = sequence::scan(Out + 2, Out + 2, n, addF<intT>(), (intT)0);
  for (long i = 0; i < n; i++) {
    intT *O = Out + (2 + n + Out[i + 2]);
    wghVertex<intT> v = G.V[i];
    for (long j = 0; j < v.degree; j++) {
      O[j] = v.Neighbors[j];
      O[j + m] = v.nghWeights[j];
    }
  }
  int r = writeArrayToFile(WghAdjGraphHeader, Out, totalLen, fname);
  free(Out);
  return r;
}

template <class intT> edgeArray<intT> readSNAP(char *fname) {
  _seq<char> S = readStringFromFile(fname);
  char *S2 = newA(char, S.n);
  // ignore starting lines with '#' and find where to start in file
  unsigned long k = 0;
  while (1) {
    if (S.A[k] == '#') {
      while (S.A[k++] != '\n')
        continue;
    }
    if (k >= S.n || S.A[k] != '#')
      break;
  }
  parallel_for(unsigned long i = 0; i < S.n - k; i++) S2[i] = S.A[k + i];
  S.del();

  words W = stringToWords(S2, S.n - k);
  unsigned long n = W.m / 2;
  edge<intT> *E = newA(edge<intT>, n);

  {
    parallel_for(unsigned long i = 0; i < n; i++) {
      // E[i] = edge<intT>(atol(W.Strings[2 * i]),
      // atol(W.Strings[2 * i + 1]));
      E[i].u = atol(W.Strings[2 * i]);
      E[i].v = atol(W.Strings[2 * i + 1]);
    }
  }
  W.del();

  unsigned long maxR = 0;
  unsigned long maxC = 0;
  for (unsigned long i = 0; i < n; i++) {
    maxR = max<intT>(maxR, E[i].u);
    maxC = max<intT>(maxC, E[i].v);
  }
  unsigned long maxrc = max<intT>(maxR, maxC) + 1;
  return edgeArray<intT>(E, maxrc, maxrc, n);
}

template <class intT> wghEdgeArray<intT> readWghSNAP(char *fname) {
  _seq<char> S = readStringFromFile(fname);
  char *S2 = newA(char, S.n);
  // ignore starting lines with '#' and find where to start in file
  long k = 0;
  while (1) {
    if (S.A[k] == '#') {
      while (S.A[k++] != '\n')
        continue;
    }
    if (k >= S.n || S.A[k] != '#')
      break;
  }
  parallel_for(long i = 0; i < S.n - k; i++) S2[i] = S.A[k + i];
  S.del();

  words W = stringToWords(S2, S.n - k);
  long n = W.m / 3;
  wghEdge<intT> *E = newA(wghEdge<intT>, n);
  {
    parallel_for(long i = 0; i < n; i++) E[i] =
        wghEdge<intT>(atol(W.Strings[3 * i]), atol(W.Strings[3 * i + 1]),
                      atol(W.Strings[3 * i + 2]));
  }
  W.del();

  long maxR = 0;
  long maxC = 0;
  for (long i = 0; i < n; i++) {
    maxR = max<intT>(maxR, E[i].u);
    maxC = max<intT>(maxC, E[i].v);
  }
  long maxrc = max<intT>(maxR, maxC) + 1;
  return wghEdgeArray<intT>(E, maxrc, maxrc, n);
}

template <class intT> graph<intT> readGraphFromFile(char *fname) {
  _seq<char> S = readStringFromFile(fname);
  words W = stringToWords(S.A, S.n);
  if (W.Strings[0] != AdjGraphHeader) {
    cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
    abort();
  }

  long len = W.m - 1;
  uintT *In = newA(uintT, len);
  { parallel_for(long i = 0; i < len; i++) In[i] = atol(W.Strings[i + 1]); }
  W.del();

  long n = In[0];
  long m = In[1];

  if (len != n + m + 2) {
    cout << "Bad input file: length = " << len << " n+m+2 = " << n + m + 2
         << endl;
    abort();
  }
  vertex<intT> *v = newA(vertex<intT>, n);
  uintT *offsets = In + 2;
  uintT *edges = In + 2 + n;

  parallel_for(uintT i = 0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = (intT *)(edges + o);
  }
  return graph<intT>(v, (intT)n, (uintT)m, (intT *)In);
}

template <class intT> wghGraph<intT> readWghGraphFromFile(char *fname) {
  _seq<char> S = readStringFromFile(fname);
  words W = stringToWords(S.A, S.n);
  if (W.Strings[0] != WghAdjGraphHeader) {
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m - 1;
  intT *In = newA(intT, len);
  { parallel_for(long i = 0; i < len; i++) In[i] = atol(W.Strings[i + 1]); }
  W.del();

  long n = In[0];
  long m = In[1];

  if (len != n + 2 * m + 2) {
    cout << "Bad input file" << endl;
    abort();
  }
  wghVertex<intT> *v = newA(wghVertex<intT>, n);
  uintT *offsets = (uintT *)In + 2;
  uintT *edges = (uintT *)In + 2 + n;
  intT *weights = In + 2 + n + m;
  parallel_for(uintT i = 0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = (intT *)(edges + o);
    v[i].nghWeights = (weights + o);
  }
  return wghGraph<intT>(v, (intT)n, (uintT)m, (intT *)In, weights);
}
}; // namespace benchIO

#endif // _BENCH_GRAPH_IO
