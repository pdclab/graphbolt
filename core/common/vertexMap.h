// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
// 
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

#ifndef __VERTEX_MAP_H__
#define __VERTEX_MAP_H__
#include "../graph/vertex.h"
#include "../graph/vertexSubset.h"
#include "../graph/IO.h"

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
const flags dense_in = 64;
const flags dense_in_pull = 128;
inline bool should_output(const flags &fl) { return !(fl & no_output); }

//*****VERTEX FUNCTIONS*****
template <class F, class VS,
          typename std::enable_if<!std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
void vertexMap(VS &V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parallel_for(long i = 0; i < n; i++) {
      if (V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }
  } else {
    parallel_for(long i = 0; i < m; i++) { f(V.vtx(i), V.vtxData(i)); }
  }
}

template <class VS, class F,
          typename std::enable_if<std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
void vertexMap(VS &V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parallel_for(long i = 0; i < n; i++) {
      if (V.isIn(i)) {
        f(i);
      }
    }
  } else {
    parallel_for(long i = 0; i < m; i++) { f(V.vtx(i)); }
  }
}

// Note: this is the version of vertexMap in which only a subset of the
// input vertexSubset is returned
template <class F> vertexSubset vertexFilter(vertexSubset &V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool *d_out = newA(bool, n);
  { parallel_for(long i = 0; i < n; i++) d_out[i] = 0; }
  { parallel_for(long i = 0; i < n; i++) if (V.d[i]) d_out[i] = filter(i); }
  // {for(long i=0;i<n;i++)
  //     if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n, d_out);
}
#endif