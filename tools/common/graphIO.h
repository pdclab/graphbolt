// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval
// Vora.
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
typedef pair<long, pair<long, char *>> intTriple;

// **************************************************************
//    EDGE ARRAY REPRESENTATION
// **************************************************************

struct edge {
  uintV u;
  uintV v;
  edge(uintV f, uintV s) : u(f), v(s) {}
  edge() : u(0), v(0) {}
};

struct edgeArray {
  edge *E;
  uintV numRows;
  uintV numCols;
  long nonZeros;
  void del() { free(E); }
  edgeArray(edge *EE, uintV r, uintV c, long nz)
      : E(EE), numRows(r), numCols(c), nonZeros(nz) {}
  edgeArray() {}
};

struct wghEdge {
  uintV u;
  uintV v;
  char *w;
  wghEdge(uintV f, uintV s, char* t) : u(f), v(s), w(t) {}
};

struct wghEdgeArray {
  wghEdge *E;
  uintV numRows;
  uintV numCols;
  long nonZeros;
  void del() { free(E); }
  wghEdgeArray(wghEdge *EE, uintV r, uintV c, long nz)
      : E(EE), numRows(r), numCols(c), nonZeros(nz) {}
  wghEdgeArray() {}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

struct vertex {
  uintV *Neighbors;
  intE degree;
  void del() {}
  vertex(uintV *N, intE d) : Neighbors(N), degree(d) {}
};

struct wghVertex {
  uintV *Neighbors;
  intE degree;
  char **nghWeights;
  void del() {
    free(Neighbors);
    free(nghWeights);
  }
  wghVertex(uintV *N, char **W, intE d)
      : Neighbors(N), nghWeights(W), degree(d) {}
};

struct graph {
  vertex *V;
  uintV n;
  long m;
  long *offsets;
  uintV *edges;

  graph(vertex *VV, uintV nn, long mm)
      : V(VV), n(nn), m(mm), offsets(NULL), edges(NULL) {}

  graph(vertex *VV, uintV nn, long mm, long *_offsets, uintV *_edges)
      : V(VV), n(nn), m(mm), offsets(_offsets), edges(_edges) {}

  void del() {
    if (offsets != NULL)
      free(offsets);
    if (edges != NULL)
      free(edges);
    free(V);
  }
};

struct wghGraph {
  wghVertex *V;
  uintV n;
  long m;
  long *offsets;
  uintV *edges;
  char **weights;

  wghGraph(wghVertex *VV, intV nn, uintV mm, long* _offsets, uintV *_edges, char **_weights)
      : V(VV), n(nn), m(mm), offsets(_offsets), edges(_edges), weights(_weights) {}

  void del() {
    if (offsets != NULL) {
      free(offsets);
    }
    if (edges != NULL) {
      free(edges);
    }
    if (weights != NULL) {
      for (int i = 0; i < m; i++) {
        free(weights[i]);
      }
      free(weights);
    }
    free(V);
  }
};

// **************************************************************
//    GRAPH UTILITIES
// **************************************************************
struct edgeCmp {
  bool operator()(edge e1, edge e2) {
    return ((e1.u < e2.u) ? 1 : ((e1.u > e2.u) ? 0 : (e1.v < e2.v)));
  }
};

struct wghEdgeCmp {
  bool operator()(wghEdge e1, wghEdge e2) {
    return ((e1.u < e2.u) ? 1 : ((e1.u > e2.u) ? 0 : (e1.v < e2.v)));
  }
};

struct edgeCopy {
  void operator()(edge e1, edge e2) {
    e2.u = e1.u;
    e2.v = e1.v;
  }
};

edgeArray remDuplicates(edgeArray A) {
  long m = A.nonZeros;
  edge *E = newA(edge, m);
  {
    parallel_for(long i = 0; i < m; i++) {
      E[i].u = A.E[i].u;
      E[i].v = A.E[i].v;
    }
  }
  quickSort(E, m, edgeCmp());
  long *flags = newA(long, m);
  flags[0] = 1;
  {
    parallel_for(long i = 1; i < m; i++) {
      if ((E[i].u != E[i - 1].u) || (E[i].v != E[i - 1].v))
        flags[i] = 1;
      else
        flags[i] = 0;
    }
  }

  long mm = sequence::plusScan(flags, flags, m);
  edge *F = newA(edge, mm);
  F[mm - 1] = E[m - 1];
  {
    parallel_for(long i = 0; i < m - 1; i++) {
      if (flags[i] != flags[i + 1])
        F[flags[i]] = E[i];
    }
  }
  free(flags);
  free(E);
  return edgeArray(F, A.numRows, A.numCols, mm);
}

wghEdgeArray remDuplicates(wghEdgeArray A) {
  long m = A.nonZeros;
  wghEdge *E = newA(wghEdge, m);
  {
    parallel_for(long i = 0; i < m; i++) {
      E[i].u = A.E[i].u;
      E[i].v = A.E[i].v;
      E[i].w = A.E[i].w;
    }
  }
  quickSort(E, m, wghEdgeCmp());
  long *flags = newA(long, m);
  flags[0] = 1;
  {
    parallel_for(long i = 1; i < m; i++) {
      if ((E[i].u != E[i - 1].u) || (E[i].v != E[i - 1].v))
        flags[i] = 1;
      else
        flags[i] = 0;
    }
  }

  long mm = sequence::plusScan(flags, flags, m);
  wghEdge *F = newA(wghEdge, mm);
  F[mm - 1] = E[m - 1];
  {
    parallel_for(long i = 0; i < m - 1; i++) {
      if (flags[i] != flags[i + 1])
        F[flags[i]] = E[i];
    }
  }
  free(flags);
  free(E);
  return wghEdgeArray(F, A.numRows, A.numCols, mm);
}

struct nEQF {
  bool operator()(edge e) { return (e.u != e.v); }
};

struct nEQFWgh {
  bool operator()(wghEdge e) { return (e.u != e.v); }
};

edgeArray makeSymmetric(edgeArray A) {
  long m = A.nonZeros;
  edge *E = A.E;
  edge *F = newA(edge, 2 * m);
  long mm = sequence::filter(E, F, m, nEQF());
  parallel_for(long i = 0; i < mm; i++) {
    F[i + mm].u = F[i].v;
    F[i + mm].v = F[i].u;
  }

  edgeArray R = remDuplicates(edgeArray(F, A.numRows, A.numCols, 2 * mm));
  free(F);

  return R;
}

wghEdgeArray makeSymmetric(wghEdgeArray A) {
  long m = A.nonZeros;
  wghEdge *E = A.E;
  wghEdge *F = newA(wghEdge, 2 * m);
  long mm = sequence::filter(E, F, m, nEQFWgh());
  parallel_for(long i = 0; i < mm; i++) {
    F[i + mm].u = F[i].v;
    F[i + mm].v = F[i].u;
    F[i + mm].w = F[i].w;
  }

  wghEdgeArray R = remDuplicates(wghEdgeArray(F, A.numRows, A.numCols, 2 * mm));
  free(F);

  return R;
}

struct getuF {
  uintV operator()(edge e) { return e.u; }
};

template <class E> struct pairFirstCmp {
  uintV operator()(E e) { return e.u; }

  bool operator()(E a, E b) { return a.u < b.u; }
};

template <class E> struct pairFirstCmp1 {
  bool operator()(pair<E, E> a, pair<E, E> b) { return a.first < b.first; }
};

graph graphFromEdges(edgeArray EA, bool makeSym) {
  edgeArray A;
  if (makeSym) {
    cout << "Make symmetric .... \n";
    A = makeSymmetric(EA);
    cout << "Make symmetric done\n";
  } else { // should have copy constructor
    edge *E = newA(edge, EA.nonZeros);
    parallel_for(long i = 0; i < EA.nonZeros; i++) E[i] = EA.E[i];
    A = edgeArray(E, EA.numRows, EA.numCols, EA.nonZeros);
  }
  long m = A.nonZeros;
  long n = max<long>(A.numCols, A.numRows);

  quickSort(A.E, m, [](edge &val1, edge &val2) { return val1.u < val2.u; });

  // long *offsets = newA(long, n * 2);
  long *offsets = newA(long, n);
  parallel_for(long i = 0; i < n; i++) { offsets[i] = m; }

  parallel_for(long i = 0; i < m - 1; i++) {
    uintV currV = A.E[i].u;
    uintV nextV = A.E[i + 1].u;
    if (currV != nextV) {
      offsets[nextV] = i + 1;
    }
  }
  offsets[A.E[0].u] = 0;
  sequence::scanIBack(offsets, offsets, (long)n, minF<long>(), (long)m);

  // #endif
  uintV *outEdges = newA(uintV, m);
  vertex *v = newA(vertex, n);

  parallel_for(uintV i = 0; i < n; i++) {
    long o = offsets[i];
    long l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = outEdges + o;
    for (long j = 0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o + j].v;
    }
    if (l > 1) {
      quickSort(v[i].Neighbors, l,
                [](uintV val1, uintV val2) { return val1 < val2; });
    }
  }
  A.del();
  EA.del();
  // free(offsets);
  return graph(v, n, m, offsets, outEdges);
}

wghGraph wghGraphFromWghEdges(wghEdgeArray EA, bool makeSym) {
  wghEdgeArray A;
  if (makeSym) {
    cout << "Make symmetric .... \n";
    A = makeSymmetric(EA);
    cout << "Make symmetric done\n";
  } else { // should have copy constructor
    wghEdge *E = newA(wghEdge, EA.nonZeros);
    parallel_for(long i = 0; i < EA.nonZeros; i++) E[i] = EA.E[i];
    A = wghEdgeArray(E, EA.numRows, EA.numCols, EA.nonZeros);
  }
  long m = A.nonZeros;
  long n = max<long>(A.numCols, A.numRows);
  
  quickSort(A.E, m, [](wghEdge &val1, wghEdge &val2) { return val1.u < val2.u; });
  
  long *offsets = newA(long, n);
  parallel_for(long i = 0; i < n; i++) { offsets[i] = m; }

  parallel_for(long i = 0; i < m - 1; i++) {
    uintV currV = A.E[i].u;
    uintV nextV = A.E[i + 1].u;
    if (currV != nextV) {
      offsets[nextV] = i + 1;
    }
  }
  offsets[A.E[0].u] = 0;
  sequence::scanIBack(offsets, offsets, (long)n, minF<long>(), (long)m);

  uintV *outEdges = newA(uintV, m);
  char **outWeights = newA(char *, m);
  wghVertex *v = newA(wghVertex, n);

  parallel_for(uintV i = 0; i < n; i++) {
    long o = offsets[i];
    long l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = outEdges + o;
    v[i].nghWeights = outWeights + o;
    for (long j = 0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o + j].v;
      v[i].nghWeights[j] = A.E[o + j].w;
    }
  }
  A.del();
  EA.del();
  // free(offsets);
  return wghGraph(v, n, m, offsets, outEdges, outWeights);
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


int writeArrayToFile(string header, long *offsets, uintV *edges, long n, long m, char *fileName) {
  
  ofstream file(fileName, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return 1;
  }
  file << header << endl;
  file << n << endl;
  file << m << endl;
  writeArrayToStream(file, offsets, n);
  writeArrayToStream(file, edges, m);
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

template <class intT> int xToStringLen(edge a) {
  return xToStringLen(a.u) + xToStringLen(a.v) + 1;
}

template <class intT> void xToString(char *s, edge a) {
  int l = xToStringLen(a.u);
  xToString(s, a.u);
  s[l] = ' ';
  xToString(s + l + 1, a.v);
}

namespace benchIO {
using namespace std;

string AdjGraphHeader = "AdjacencyGraph";
string WghAdjGraphHeader = "WeightedAdjacencyGraph";

int writeGraphToFile(graph G, char *fname) {
  long m = G.m;
  long n = G.n;
  long totalLen = 2 + n + m;
  long *offsets = newA(long, n);
  uintV *edges = newA(uintV, m);

  parallel_for(long i = 0; i < n; i++) { offsets[i] = G.V[i].degree; }
  long total = sequence::scan(offsets, offsets, n, addF<long>(), (long)0);

  for (long i = 0; i < n; i++) {
    uintV *O = edges + offsets[i];

    vertex v = G.V[i];
    for (intE j = 0; j < v.degree; j++)
      O[j] = v.Neighbors[j];
  }

  int r = writeArrayToFile(AdjGraphHeader, offsets, edges, n, m, fname);
  free(offsets);
  free(edges);
  return r;
}

int writeWghGraphToFile(wghGraph G, char *fname) {
  ofstream file(fname, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fname << std::endl;
    return 1;
  }
  file << WghAdjGraphHeader << endl;
  file << G.n << endl;
  file << G.m << endl;
  writeArrayToStream(file, G.offsets, G.n);
  writeArrayToStream(file, G.edges, G.m);
  writeArrayToStream(file, G.weights, G.m);
  file.close();
  return 0;
}

edgeArray readSNAP(char *fname) {
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
  edge *E = newA(edge, n);

  {
    parallel_for(unsigned long i = 0; i < n; i++) {
      // E[i] = edge(atol(W.Strings[2 * i]),
      // atol(W.Strings[2 * i + 1]));
      E[i].u = atol(W.Strings[2 * i]);
      E[i].v = atol(W.Strings[2 * i + 1]);
    }
  }
  W.del();

  unsigned long maxR = 0;
  unsigned long maxC = 0;
  for (unsigned long i = 0; i < n; i++) {
    maxR = max<uintV>(maxR, E[i].u);
    maxC = max<uintV>(maxC, E[i].v);
  }
  unsigned long maxrc = max<uintV>(maxR, maxC) + 1;
  return edgeArray(E, maxrc, maxrc, n);
}

wghEdgeArray readWghSNAP(char *fname) {
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
  unsigned long n = W.m / 3;
  wghEdge *E = newA(wghEdge, n);
  {
    parallel_for(unsigned long i = 0; i < n; i++) {
      E[i] = wghEdge(atol(W.Strings[3 * i]), atol(W.Strings[3 * i + 1]),
                     W.Strings[3 * i + 2]);
    }
  }

  unsigned long maxR = 0;
  unsigned long maxC = 0;
  for (long i = 0; i < n; i++) {
    maxR = max<uintV>(maxR, E[i].u);
    maxC = max<uintV>(maxC, E[i].v);
  }
  unsigned long maxrc = max<uintV>(maxR, maxC) + 1;
  return wghEdgeArray(E, maxrc, maxrc, n);
}

graph readGraphFromFile(char *fname) {
  _seq<char> S = readStringFromFile(fname);
  words W = stringToWords(S.A, S.n);
  if (W.Strings[0] != AdjGraphHeader) {
    cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
    abort();
  }

  long len = W.m - 1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);

  long *offsets = newA(long, n);
  uintV *outEdges = newA(uintV, m);
  { parallel_for(long i = 0; i < n; i++) offsets[i] = atol(W.Strings[i + 2 + 1]); }
  { parallel_for(long i = 0; i < m; i++) outEdges[i] = atol(W.Strings[i + 2 + n + 1]); }

  W.del();

  if (len != n + m + 2) {
    cout << "Bad input file: length = " << len << " n+m+2 = " << n + m + 2
         << endl;
    abort();
  }

  vertex *v = newA(vertex, n);

  parallel_for(long i = 0; i < n; i++) {
    long o = offsets[i];
    long l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].degree = l;
    v[i].Neighbors = (outEdges + o);
  }
  return graph(v, n, m, offsets, outEdges);
}


}; // namespace benchIO

#endif // _BENCH_GRAPH_IO
