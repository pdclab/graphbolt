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
  uintV *neighbors;
#ifdef EDGEDATA
  EdgeData *edgeDataArray;
#endif
  intE degree;
  void del() { free(neighbors); }
  symmetricVertex(uintV *n, intE d) : neighbors(n), degree(d) {}
  uintV *getInNeighbors() { return neighbors; }
  const uintV *getInNeighbors() const { return neighbors; }
  uintV *getOutNeighbors() { return neighbors; }
  const uintV *getOutNeighbors() const { return neighbors; }
  uintV getInNeighbor(intE j) const { return neighbors[j]; }
  uintV getOutNeighbor(intE j) const { return neighbors[j]; }

  void setInNeighbor(intE j, uintV ngh) { neighbors[j] = ngh; }
  void setOutNeighbor(intE j, uintV ngh) { neighbors[j] = ngh; }
  void setInNeighbors(uintV *_i) { neighbors = _i; }
  void setOutNeighbors(uintV *_i) { neighbors = _i; }

#ifdef EDGEDATA
  EdgeData *getInEdgeData(intE j) const { return &edgeDataArray[j]; }
  EdgeData *getOutEdgeData(intE j) const { return &edgeDataArray[j]; }
  EdgeData *getInEdgeDataArray() const { return edgeDataArray; }
  EdgeData *getOutEdgeDataArray() const { return edgeDataArray; }
  void setInEdgeDataElement(intE j, EdgeData data) { edgeDataArray[j] = data; }
  void setOutEdgeDataElement(intE j, EdgeData data) {
    edgeDataArray[j] = data;
  }
  void setInEdgeDataArray(EdgeData *_i) { edgeDataArray = _i; }
  void setOutEdgeDataArray(EdgeData *_i) { edgeDataArray = _i; }
#endif

  intE getInDegree() const { return degree; }
  intE getOutDegree() const { return degree; }
  void setInDegree(intE _d) { degree = _d; }
  void setOutDegree(intE _d) { degree = _d; }
  intE fetchAndAddInDegree(int delta) {
    intE toReturn = pbbs::fetch_and_add(&degree, delta);
    return toReturn;
  }
  intE fetchAndAddOutDegree(int delta) {
    intE toReturn = pbbs::fetch_and_add(&degree, delta);
    return toReturn;
  }
};

struct asymmetricVertex {
  uintV *inNeighbors, *outNeighbors;
#ifdef EDGEDATA
  EdgeData *outEdgeDataArray, *inEdgeDataArray;
#endif
  intE outDegree;
  intE inDegree;
  void del() {
    free(inNeighbors);
    free(outNeighbors);
  }
  asymmetricVertex(uintV *iN, uintV *oN, intE id, intE od)
      : inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {}
  uintV *getInNeighbors() { return inNeighbors; }
  const uintV *getInNeighbors() const { return inNeighbors; }
  uintV *getOutNeighbors() { return outNeighbors; }
  const uintV *getOutNeighbors() const { return outNeighbors; }
  uintV getInNeighbor(intE j) const { return inNeighbors[j]; }
  uintV getOutNeighbor(intE j) const { return outNeighbors[j]; }
  void setInNeighbor(intE j, uintV ngh) { inNeighbors[j] = ngh; }
  void setOutNeighbor(intE j, uintV ngh) { outNeighbors[j] = ngh; }
  void setInNeighbors(uintV *_i) { inNeighbors = _i; }
  void setOutNeighbors(uintV *_i) { outNeighbors = _i; }
  intE getInDegree() const { return inDegree; }
  intE getOutDegree() const { return outDegree; }
  void setInDegree(intE _d) { inDegree = _d; }
  void setOutDegree(intE _d) { outDegree = _d; }

#ifdef EDGEDATA
  EdgeData *getInEdgeData(intE j) const { return &inEdgeDataArray[j]; }
  EdgeData *getOutEdgeData(intE j) const { return &outEdgeDataArray[j]; }
  EdgeData *getInEdgeDataArray() const { return inEdgeDataArray; }
  EdgeData *getOutEdgeDataArray() const { return outEdgeDataArray; }
  void setInEdgeDataElement(intE j, EdgeData data) {
    inEdgeDataArray[j] = data;
  }
  void setOutEdgeDataElement(intE j, EdgeData data) {
    outEdgeDataArray[j] = data;
  }
  void setInEdgeDataArray(EdgeData *_i) { inEdgeDataArray = _i; }
  void setOutEdgeDataArray(EdgeData *_i) { outEdgeDataArray = _i; }
#endif

  intE fetchAndAddInDegree(int delta) {
    intE toReturn = pbbs::fetch_and_add(&inDegree, delta);
    return toReturn;
  }
  intE fetchAndAddOutDegree(int delta) {
    intE toReturn = pbbs::fetch_and_add(&outDegree, delta);
    return toReturn;
  }
};

#endif
