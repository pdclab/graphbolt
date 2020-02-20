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

#ifndef VERTEX_H
#define VERTEX_H
#include "vertexSubset.h"
using namespace std;


struct symmetricVertex {
  uintE *neighbors;
  uintT degree;
  void del() { free(neighbors); }
  symmetricVertex(uintE *n, uintT d)
      : neighbors(n), degree(d) {
  }
  uintE *getInNeighbors() { return neighbors; }
  const uintE *getInNeighbors() const { return neighbors; }
  uintE *getOutNeighbors() { return neighbors; }
  const uintE *getOutNeighbors() const { return neighbors; }
  uintE getInNeighbor(uintT j) const { return neighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return neighbors[j]; }

  void setInNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
  void setInNeighbors(uintE *_i) { neighbors = _i; }
  void setOutNeighbors(uintE *_i) { neighbors = _i; }

  uintT getInDegree() const { return degree; }
  uintT getOutDegree() const { return degree; }
  void setInDegree(uintT _d) { degree = _d; }
  void setOutDegree(uintT _d) { degree = _d; }
  uintT fetchAndAddInDegree(int delta) {
    uintT toReturn = pbbs::fetch_and_add(&degree, delta);
    return toReturn;
  }
  uintT fetchAndAddOutDegree(int delta) {
    uintT toReturn = pbbs::fetch_and_add(&degree, delta);
    return toReturn;
  }

};

struct asymmetricVertex {
  uintE *inNeighbors, *outNeighbors;
  uintT outDegree;
  uintT inDegree;
  void del() {
    free(inNeighbors);
    free(outNeighbors);
  }
  asymmetricVertex(uintE *iN, uintE *oN, uintT id, uintT od)
      : inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {
  }
  uintE *getInNeighbors() { return inNeighbors; }
  const uintE *getInNeighbors() const { return inNeighbors; }
  uintE *getOutNeighbors() { return outNeighbors; }
  const uintE *getOutNeighbors() const { return outNeighbors; }
  uintE getInNeighbor(uintT j) const { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setInNeighbors(uintE *_i) { inNeighbors = _i; }
  void setOutNeighbors(uintE *_i) { outNeighbors = _i; }
  uintT getInDegree() const { return inDegree; }
  uintT getOutDegree() const { return outDegree; }
  void setInDegree(uintT _d) { inDegree = _d; }
  void setOutDegree(uintT _d) { outDegree = _d; }
  uintT fetchAndAddInDegree(int delta) {
    uintT toReturn = pbbs::fetch_and_add(&inDegree, delta);
    return toReturn;
  }
  uintT fetchAndAddOutDegree(int delta) {
    uintT toReturn = pbbs::fetch_and_add(&outDegree, delta);
    return toReturn;
  }

};

#endif
