// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
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

#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H
#include "graph.h"
#include "vertexSubset.h"

template <class vertex>
uintE countOutNghs(vertex *V, vertexSubset vs, uintE n) {
  uintE count = 0;
  for (uintE i = 0; i < n; i++) {
    if (vs.isIn(i)) {
      // cout << i << count << ", ";
      count += V[i].getOutDegree();
      // cout << count << "\n";
    }
  }
  return count;
}

template <class vertex, class check>
uintE countOutNghs(vertex *V, vertexSubset vs, check checkf, uintE n) {
  uintE count = 0;
  for (uintE i = 0; i < n; i++) {
    if (vs.isIn(i)) {
      for (uintE j = 0; j < V[i].getOutDegree(); j++) {
        // Check if the outNgh is valid before incrementing the count
        uintE outNgh = V[i].getOutNeighbor(j);
        if (checkf(outNgh)) {
          count++;
        }
      }
    }
  }
  return count;
}

#endif
