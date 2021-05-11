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

#ifndef GRAPH_H
#define GRAPH_H
#include "../common/parallel.h"
#include "../common/quickSort.h"
#include "vertex.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <unordered_set>
#include <vector>

using namespace std;

struct edge {
  uintV source;
  uintV destination;
#ifdef EDGEDATA
  EdgeData *edgeData; // pointer to struct
  edge(uintV f, uintV s, EdgeData *e) 
      : source(f), destination(s), edgeData(e) {}
#else
  edge(uintV f, uintV s) : source(f), destination(s) {}
#endif
};

struct edgeArray {
  edge *E;
#ifdef EDGEDATA
  EdgeData *edgeDataArray; // (pointer to array of structs)
#endif
  long size;
  long maxVertex;
  void del() {
    if (E != nullptr) {
#ifdef EDGEDATA
      parallel_for(uintV i = 0; i < size; i++) { edgeDataArray[i].del(); }
      free(edgeDataArray);
#endif
      free(E);
    }
  }
  edgeArray() : size(0) { E = nullptr; }
#ifdef EDGEDATA
  edgeArray(edge *EE, EdgeData *_edgeDataArray, long _size, long _maxVertex)
      : E(EE), edgeDataArray(_edgeDataArray), size(_size),
        maxVertex(_maxVertex) {}
#else
  edgeArray(edge *EE, long _size, long _maxVertex)
      : E(EE), size(_size), maxVertex(_maxVertex) {}
#endif
};

struct edgesToDelete {
  vector<uintV> outEdgesToDelete;
  vector<uintV> inEdgesToDelete;
#ifdef EDGEDATA
  vector<EdgeData *> outEdgeDataToDelete;
  vector<EdgeData *> inEdgeDataToDelete;
#endif

  edgesToDelete() {}
  edgesToDelete(vector<uintV> _outEdgesToDelete, vector<uintV> _inEdgesToDelete)
      : outEdgesToDelete(_outEdgesToDelete), inEdgesToDelete(_inEdgesToDelete) {
  }

  inline bool hasOutNeighbor(uintV outNeighbor) {
    vector<uintV>::iterator index =
        find(outEdgesToDelete.begin(), outEdgesToDelete.end(), outNeighbor);
    if (index != outEdgesToDelete.end()) {
      return true;
    }
    return false;
  }

  inline bool hasInNeighbor(uintV inNeighbor) {
    vector<uintV>::iterator index =
        find(inEdgesToDelete.begin(), inEdgesToDelete.end(), inNeighbor);
    if (index != inEdgesToDelete.end()) {
      return true;
    }
    return false;
  }

  inline long getIndexOfOutNeighbor(uintV outNeighbor, long start = 0) {

    vector<uintV>::iterator it = find(outEdgesToDelete.begin() + start,
                                      outEdgesToDelete.end(), outNeighbor);
    if (it != outEdgesToDelete.end()) {
      return it - outEdgesToDelete.begin();
    }
    return -1;
  }

  inline long getIndexOfInNeighbor(uintV inNeighbor, long start = 0) {
    vector<uintV>::iterator it = find(inEdgesToDelete.begin() + start,
                                      inEdgesToDelete.end(), inNeighbor);
    if (it != inEdgesToDelete.end()) {
      return it - inEdgesToDelete.begin();
    }
    return -1;
  }

  inline void clear() {
    outEdgesToDelete.clear();
    inEdgesToDelete.clear();
#ifdef EDGEDATA
    outEdgeDataToDelete.clear();
    inEdgeDataToDelete.clear();
#endif
  }

  inline void insertOutEdge(uintV v) { outEdgesToDelete.push_back(v); }

  inline void insertInEdge(uintV v) { inEdgesToDelete.push_back(v); }

#ifdef EDGEDATA
  inline void insertOutData(EdgeData *edgeData) {
    outEdgeDataToDelete.push_back(edgeData);
  }

  inline void insertInData(EdgeData *edgeData) {
    inEdgeDataToDelete.push_back(edgeData);
  }
#endif
};

template <class E> struct SimpleCmp {
  bool operator()(E a, E b) { return a < b; }
};

template <class edge> struct MinSrcCmp {
  bool operator()(edge a, edge b) { return a.source < b.source; }
};

template <class edge> struct MinDesCmp {
  bool operator()(edge a, edge b) { return a.destination < b.destination; }
};

struct edgeDeletionData {
  edgesToDelete *dataMap;
  unsigned long numberOfDeletions;
  edge *edgesArray;
  uintV n;
  bool *updatedVertices;

  edgeDeletionData(uintV _n) : n(_n) {
    updatedVertices = newA(bool, n);
    dataMap = new edgesToDelete[n];
    parallel_for(uintV i = 0; i < n; i++) { updatedVertices[i] = 0; }
    numberOfDeletions = 0;
    edgesArray = nullptr;
  }

  inline void addEdge(edge &E) {
    dataMap[E.source].insertOutEdge(E.destination);
    dataMap[E.destination].insertInEdge(E.source);
    updatedVertices[E.source] = 1;
    updatedVertices[E.destination] = 1;
  }

  void updateWithEdgesArray(edgeArray &edgeArrayToDelete) {
    edgesArray = edgeArrayToDelete.E;
    numberOfDeletions = edgeArrayToDelete.size;

    // sort on src
    quickSort(edgesArray, numberOfDeletions, MinSrcCmp<edge>());

    parallel_for(long i = 0; i < numberOfDeletions; i++) {
      long prev = (i == 0) ? 0 : i - 1;
      if ((i == 0) || (edgesArray[prev].source != edgesArray[i].source)) {
        long j = i;
        updatedVertices[edgesArray[i].source] = 1;
        while ((j < numberOfDeletions) &&
               (edgesArray[i].source == edgesArray[j].source)) {
          dataMap[edgesArray[j].source].insertOutEdge(
              edgesArray[j].destination);
#ifdef EDGEDATA
          dataMap[edgesArray[j].source].insertOutData(edgesArray[j].edgeData);
#endif
          j++;
        }
      }
    }

    // sort on dest
    quickSort(edgesArray, numberOfDeletions, MinDesCmp<edge>());

    parallel_for(long i = 0; i < numberOfDeletions; i++) {
      long prev = (i == 0) ? 0 : i - 1;
      if (i == 0 || edgesArray[prev].destination != edgesArray[i].destination) {
        long j = i;
        updatedVertices[edgesArray[i].destination] = 1;
        while (j < numberOfDeletions &&
               edgesArray[i].destination == edgesArray[j].destination) {
          dataMap[edgesArray[j].destination].insertInEdge(edgesArray[j].source);
#ifdef EDGEDATA
          dataMap[edgesArray[j].destination].insertInData(
              edgesArray[j].edgeData);
#endif
          j++;
        }
      }
    }
  }

  edgesToDelete &getEdgeDeletionData(uintV vertex) { return dataMap[vertex]; }

  void reset() {
    numberOfDeletions = 0;
    parallel_for(uintV i = 0; i < n; i++) {
      if (updatedVertices[i] == 1) {
        // Delete the entries in dataMap corresponding to the vertex i
        dataMap[i].clear();
        updatedVertices[i] = 0;
      }
    }
  }

  void updateNumVertices(uintV maxVertex) {
    if (n) {
      delete[] dataMap;
    }
    updatedVertices = renewA(bool, updatedVertices, maxVertex);
    parallel_for(uintV i = 0; i < maxVertex; i++) { updatedVertices[i] = 0; }
    dataMap = new edgesToDelete[maxVertex];
    n = maxVertex;
  }

  void del() {
    if (n) {
      free(updatedVertices);
      delete[] dataMap;
    }
  }
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

// Class that handles implementation specific freeing of memory
// owned by the graph
struct Deletable {
public:
  virtual uintEE getNumberOfEdges() = 0;
  virtual void *getOutEdges() = 0;
  virtual void *getInEdges() = 0;
  virtual uintE *getOutEdgeOffsets() = 0;
  virtual uintE *getInEdgeOffsets() = 0;
  virtual void del() = 0;
  virtual void *updateVertices(uintV _verticesSize) = 0;
  virtual edgeArray addEdges(edgeArray &edgesToAdd, bool *updatedVertices) = 0;
  virtual void setSymmetric(bool flag) = 0;
  virtual edgeArray deleteEdges(edgeDeletionData &data, bool *updatedVertices,
                                bool debugFlag) = 0;
  virtual ~Deletable() = default;
};

template <class vertex> struct AdjacencyRep : public Deletable {
public:
  vertex *V;
  unsigned long n;
  unsigned long m;
  uintV **outEdges = NULL, **inEdges = NULL;
  uintE *outEdgesArraySize = NULL, *inEdgesArraySize = NULL;
#ifdef EDGEDATA
  EdgeData **outEdgeData = NULL, **inEdgeData = NULL;
#endif
  // Indicates whether the graph is directed or undirected. TRUE -> indicates
  // undirected graph
  bool symmetric;

  uintV *inEdgeUpdates = NULL;
  uintV *outEdgeUpdates = NULL;

#ifdef EDGEDATA
  AdjacencyRep(vertex *VV, unsigned long nn, unsigned long mm, uintV *ai,
               uintV *_inEdges = NULL, intE *_outEdgeOffsets = NULL,
               intE *_inEdgeOffsets = NULL, EdgeData *_outEdgeData = NULL,
               EdgeData *_inEdgeData = NULL)
#else
  AdjacencyRep(vertex *VV, unsigned long nn, unsigned long mm, uintV *ai,
               uintV *_inEdges = NULL, intE *_outEdgeOffsets = NULL,
               intE *_inEdgeOffsets = NULL)
#endif
      : V(VV), n(nn), m(mm), symmetric(false) {
    outEdges = newA(uintV *, nn);
    outEdgesArraySize = newA(uintE, nn);
    outEdgeUpdates = newA(uintV, nn);

#ifdef EDGEDATA
    if (_outEdgeData != NULL) {
      outEdgeData = newA(EdgeData *, nn);
    }
#endif

    if (_inEdges != NULL) {
      inEdges = newA(uintV *, nn);
      inEdgesArraySize = newA(uintE, nn);
      inEdgeUpdates = newA(uintV, nn);
#ifdef EDGEDATA
      inEdgeData = newA(EdgeData *, nn);
#endif
    }

    parallel_for(uintV i = 0; i < nn; i++) {
      uintE outEdgesSize = V[i].getOutDegree();
      outEdges[i] = newA(uintV, outEdgesSize);
      outEdgesArraySize[i] = outEdgesSize;

#ifdef EDGEDATA
      outEdgeData[i] = newA(EdgeData, outEdgesSize);
#endif
      for (uintE j = 0; j < outEdgesSize; j++) {
        outEdges[i][j] = ai[_outEdgeOffsets[i] + j];
#ifdef EDGEDATA
        new (outEdgeData[i] + j) EdgeData();
        outEdgeData[i][j].setEdgeDataFromPtr(
            &_outEdgeData[_outEdgeOffsets[i] + j]);
#endif
      }
      V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
      V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif

      if (_inEdges != NULL) {
        uintE inEdgesSize = V[i].getInDegree();
        inEdges[i] = newA(uintV, inEdgesSize);
        inEdgesArraySize[i] = inEdgesSize;

#ifdef EDGEDATA
        inEdgeData[i] = newA(EdgeData, inEdgesSize);
#endif
        for (uintE j = 0; j < inEdgesSize; j++) {
          inEdges[i][j] = _inEdges[_inEdgeOffsets[i] + j];
#ifdef EDGEDATA
          new (inEdgeData[i] + j) EdgeData();
          inEdgeData[i][j].setEdgeDataFromPtr(
              &_inEdgeData[_inEdgeOffsets[i] + j]);
#endif
        }
        V[i].setInNeighbors(inEdges[i]);
#ifdef EDGEDATA
        V[i].setInEdgeDataArray(inEdgeData[i]);
#endif
      }
    }

    if (ai != NULL) {
      free(ai);
      free(_outEdgeOffsets);
#ifdef EDGEDATA
      parallel_for(uintE i = 0; i < mm; i++) { _outEdgeData[i].del(); }
      free(_outEdgeData);
#endif
    }

    if (_inEdges != NULL) {
      free(_inEdges);
      free(_inEdgeOffsets);
#ifdef EDGEDATA
      parallel_for(uintE i = 0; i < mm; i++) { _inEdgeData[i].del(); }
      free(_inEdgeData);
#endif
    }
  }

  void setSymmetric(bool flag) { symmetric = flag; }

  void *getOutEdges() { return outEdges; }

  void *getInEdges() { return inEdges; }

  uintE *getOutEdgeOffsets() {
    std::cout << "Not defined" << std::endl;
    return nullptr;
  }

  uintE *getInEdgeOffsets() {
    std::cout << "Not defined" << std::endl;
    return nullptr;
  }

  bool isSymmetric() { return symmetric; }

  void *updateVertices(uintV maxVertex) {
    if (isSymmetric()) {
      return updateVertices_symmetric(maxVertex);
    }

    if (maxVertex >= n) {
      uintV currentVertexSize = n;
      n = maxVertex + 1;
      V = renewA(vertex, V, n);
      outEdges = renewA(uintV *, outEdges, n);
      inEdges = renewA(uintV *, inEdges, n);
      outEdgesArraySize = renewA(uintE, outEdgesArraySize, n);
      inEdgesArraySize = renewA(uintE, inEdgesArraySize, n);
      outEdgeUpdates = renewA(uintV, outEdgeUpdates, n);
      inEdgeUpdates = renewA(uintV, inEdgeUpdates, n);
#ifdef EDGEDATA
      outEdgeData = renewA(EdgeData *, outEdgeData, n);
      inEdgeData = renewA(EdgeData *, inEdgeData, n);
#endif

      parallel_for(uintV i = 0; i < currentVertexSize; i++) {
        V[i].setOutNeighbors(outEdges[i]);
        V[i].setInNeighbors(inEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
        V[i].setInEdgeDataArray(inEdgeData[i]);
#endif
      }

      parallel_for(uintV i = currentVertexSize; i < n; i++) {
        V[i].setOutDegree(0);
        V[i].setInDegree(0);
        // TODO : What is this doing??
        outEdges[i] = newA(uintV, 0);
        inEdges[i] = newA(uintV, 0);
        V[i].setOutNeighbors(outEdges[i]);
        V[i].setInNeighbors(inEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
        V[i].setInEdgeDataArray(inEdgeData[i]);
#endif
        outEdgesArraySize[i] = 0;
        inEdgesArraySize[i] = 0;
      }
      return V;
    }
  }

  void *updateVertices_symmetric(uintV maxVertex) {
    if (maxVertex > n) {
      uintV currentVertexSize = n;
      n = maxVertex + 1;
      V = renewA(vertex, V, n);
      outEdges = renewA(uintV *, outEdges, n);
      outEdgesArraySize = renewA(uintE, outEdgesArraySize, n);
      outEdgeUpdates = renewA(uintV, outEdgeUpdates, n);
#ifdef EDGEDATA
      outEdgeData = renewA(EdgeData *, outEdgeData, n);
#endif

      parallel_for(uintV i = 0; i < currentVertexSize; i++) {
        V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif
      }

      parallel_for(uintV i = currentVertexSize; i < n; i++) {
        V[i].setOutDegree(0);
        outEdges[i] = newA(uintV, 0);
        V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif
        outEdgesArraySize[i] = 0;
      }
      return V;
    }
  }

  edgeArray addEdges_symmetric(edgeArray &edgesToAdd, bool *updatedVertices) {
    parallel_for(uintV i = 0; i < n; i++) { outEdgeUpdates[i] = 0; }

    edge *E = edgesToAdd.E;
    uintE size = edgesToAdd.size;
    // for each new edge, we increment the outDegreeOffset of the source
    // and we increment the inDegreeOffset of the destination
    parallel_for(uintE i = 0; i < size; i++) {
      uintV source = E[i].source;
      uintV destination = E[i].destination;

      pbbs::write_add(&outEdgeUpdates[source], 1);

      pbbs::write_add(&outEdgeUpdates[destination], 1);

      updatedVertices[source] = 1;
      updatedVertices[destination] = 1;
    }

    parallel_for(uintE i = 0; i < n; i++) {
#ifdef INCLUDEEXTRABUFFERSPACE
      if (outEdgeUpdates[i] &&
          (outEdgeUpdates[i] + V[i].getOutDegree()) > outEdgesArraySize[i]) {
        outEdgesArraySize[i] =
            V[i].getOutDegree() + outEdgeUpdates[i] + EDGE_ARRAY_BUFFER;
        outEdges[i] = renewA(uintV, outEdges[i], outEdgesArraySize[i]);
        V[i].setOutNeighbors(outEdges[i]);
      }
#else
      if (outEdgeUpdates[i]) {
        outEdges[i] = renewA(uintV, outEdges[i],
                             (V[i].getOutDegree() + outEdgeUpdates[i]));
        V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
        outEdgeData[i] = renewA(EdgeData, outEdgeData[i],
                                (V[i].getOutDegree() + outEdgeUpdates[i]));
        V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif
      }
#endif
    }

    parallel_for(uintE i = 0; i < size; i++) {
      uintV source = E[i].source;
      uintV destination = E[i].destination;
#ifdef EDGEDATA
      EdgeData *edgeData = E[i].edgeData;
#endif
      vertex &sourceVertex = V[source];
      vertex &destinationVertex = V[destination];
      long outIndex = sourceVertex.fetchAndAddOutDegree(1);
      outEdges[source][outIndex] = destination;
      long inIndex = destinationVertex.fetchAndAddOutDegree(1);
      outEdges[destination][inIndex] = source;
#ifdef EDGEDATA
      new (outEdgeData[source] + outIndex) EdgeData();
      outEdgeData[source][outIndex].setEdgeDataFromPtr(edgeData);
      new (outEdgeData[destination] + inIndex) EdgeData();
      outEdgeData[destination][inIndex].setEdgeDataFromPtr(edgeData);
#endif
    }

    long curr_size = edgesToAdd.size;
#ifdef EDGEDATA
    EdgeData *edgeDataWeight = newA(EdgeData, curr_size * 2);
#endif
    edge *symE = newA(edge, curr_size * 2);
    long j = 0;
    for (long i = 0; i < curr_size; i++) {
      symE[j].source = edgesToAdd.E[i].source;
      symE[j].destination = edgesToAdd.E[i].destination;
#ifdef EDGEDATA
      new (edgeDataWeight + j) EdgeData();
      symE[j].edgeData = &edgeDataWeight[j];
      symE[j].edgeData->setEdgeDataFromPtr(&edgesToAdd.edgeDataArray[i]);
#endif
      j++;
      symE[j].source = edgesToAdd.E[i].destination;
      symE[j].destination = edgesToAdd.E[i].source;
#ifdef EDGEDATA
      new (edgeDataWeight + j) EdgeData();
      symE[j].edgeData = &edgeDataWeight[j];
      symE[j].edgeData->setEdgeDataFromPtr(&edgesToAdd.edgeDataArray[i]);
#endif
      j++;
    }
    uintV max_vertex = edgesToAdd.maxVertex;
    edgesToAdd.del();
#ifdef EDGEDATA
    edgesToAdd = edgeArray(symE, edgeDataWeight, j, max_vertex);
#else
    edgesToAdd = edgeArray(symE, j, max_vertex);
#endif

    uintE newSize = m + (edgesToAdd.size);
    m = newSize;
    return edgesToAdd;
  }

  edgeArray addEdges(edgeArray &edgesToAdd, bool *updatedVertices) {
    if (isSymmetric()) {
      return addEdges_symmetric(edgesToAdd, updatedVertices);
    }

    parallel_for(uintV i = 0; i < n; i++) {
      outEdgeUpdates[i] = 0;
      inEdgeUpdates[i] = 0;
    }

    edge *E = edgesToAdd.E;
    uintE size = edgesToAdd.size;
    // for each new edge, we increment the outDegreeOffset of the source
    // and we increment the inDegreeOffset of the destination
    parallel_for(uintE i = 0; i < size; i++) {
      uintV source = E[i].source;
      uintV destination = E[i].destination;

      pbbs::write_add(&outEdgeUpdates[source], 1);

      pbbs::write_add(&inEdgeUpdates[destination], 1);

      updatedVertices[source] = 1;
      updatedVertices[destination] = 1;
    }

    parallel_for(uintV i = 0; i < n; i++) {
#ifdef INCLUDEEXTRABUFFERSPACE
      if (outEdgeUpdates[i] &&
          (outEdgeUpdates[i] + V[i].getOutDegree()) > outEdgesArraySize[i]) {
        outEdgesArraySize[i] =
            V[i].getOutDegree() + outEdgeUpdates[i] + EDGE_ARRAY_BUFFER;
        outEdges[i] = renewA(uintV, outEdges[i], outEdgesArraySize[i]);
        V[i].setOutNeighbors(outEdges[i]);
      }
      if (inEdgeUpdates[i] &&
          (inEdgeUpdates[i] + V[i].getInDegree()) > inEdgesArraySize[i]) {
        inEdgesArraySize[i] =
            V[i].getInDegree() + inEdgeUpdates[i] + EDGE_ARRAY_BUFFER;
        inEdges[i] = renewA(uintV, inEdges[i], inEdgesArraySize[i]);
        V[i].setInNeighbors(inEdges[i]);
      }
#else
      if (outEdgeUpdates[i]) {
        outEdges[i] = renewA(uintV, outEdges[i],
                             (V[i].getOutDegree() + outEdgeUpdates[i]));
        V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
        outEdgeData[i] = renewA(EdgeData, outEdgeData[i],
                                (V[i].getOutDegree() + outEdgeUpdates[i]));
        V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif
      }
      if (inEdgeUpdates[i]) {
        inEdges[i] =
            renewA(uintV, inEdges[i], (V[i].getInDegree() + inEdgeUpdates[i]));
        V[i].setInNeighbors(inEdges[i]);
#ifdef EDGEDATA
        inEdgeData[i] = renewA(EdgeData, inEdgeData[i],
                               (V[i].getInDegree() + inEdgeUpdates[i]));
        V[i].setInEdgeDataArray(inEdgeData[i]);
#endif
      }
#endif
    }

    parallel_for(uintE i = 0; i < edgesToAdd.size; i++) {
      uintV source = E[i].source;
      uintV destination = E[i].destination;
#ifdef EDGEDATA
      EdgeData *edgeData = E[i].edgeData;
#endif
      vertex &sourceVertex = V[source];
      vertex &destinationVertex = V[destination];

      long outIndex = sourceVertex.fetchAndAddOutDegree(1);
      outEdges[source][outIndex] = destination;
      long inIndex = destinationVertex.fetchAndAddInDegree(1);
      inEdges[destination][inIndex] = source;
#ifdef EDGEDATA
      new (outEdgeData[source] + outIndex) EdgeData();
      outEdgeData[source][outIndex].setEdgeDataFromPtr(edgeData);
      new (inEdgeData[destination] + inIndex) EdgeData();
      inEdgeData[destination][inIndex].setEdgeDataFromPtr(edgeData);
#endif
    }

    unsigned long newSize = m + edgesToAdd.size;
    m = newSize;
    return edgesToAdd;
  }

  // Todo: change type of edges to be signed
  edgeArray deleteEdges_symmetric(edgeDeletionData &deletionsData,
                                  bool *updatedVertices, bool debugFlag) {
    uintV maxValue = numeric_limits<uintV>::max();
    intE numberOfSuccessfulDeletions = 0;

    edge *ED = newA(edge, deletionsData.numberOfDeletions * 2);
#ifdef EDGEDATA
    EdgeData *edgeDataWeight =
        newA(EdgeData, deletionsData.numberOfDeletions * 2);
#endif

    // TODO : Why?
    uintE *outDegree = newA(uintE, n);
    uintE *inDegree = newA(uintE, n);

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        outDegree[i] = V[i].getOutDegree();
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;

        parallel_for(intE j = 0; j < outEdgesToDelete.size(); j++) {
          uintV targetOutNgh = outEdgesToDelete[j];
          uintV *currOutEdges = outEdges[i];

          bool deletionSuccessful = false;
          for (intE k = 0; k < outDegree[i]; k++) {
            if (targetOutNgh == currOutEdges[k]) {
              bool casSuccessful;
              do {
                casSuccessful = CAS(&currOutEdges[k], targetOutNgh, maxValue);
              } while (currOutEdges[k] != maxValue);
              if (casSuccessful) {
                deletionSuccessful = true;
                break;
              }
            }
          }
          if (deletionSuccessful == false) {
            outEdgesToDelete[j] = maxValue;
          }
        }
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        uintV *currOutEdges = outEdges[i];
#ifdef EDGEDATA
        EdgeData *currOutEdgeData = outEdgeData[i];
#endif

        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        intE last_non_deleted_index;
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;
        intE to_delete_count = outEdgesToDelete.size();

        intE actual_to_delete_count = 0;
        for (intE i = 0; i < to_delete_count; i++) {
          if (outEdgesToDelete[i] != maxValue) {
            actual_to_delete_count++;
          }
        }
        intE total_swapped = 0;

        for (intE k = outDegree[i] - 1,
                  last_non_deleted_index = outDegree[i] - 1;
             k >= 0; k--) {
          if (currOutEdges[k] == maxValue) {
            currOutEdges[k] = currOutEdges[last_non_deleted_index];
#ifdef EDGEDATA
            currOutEdgeData[k].del();
            currOutEdgeData[k] = currOutEdgeData[last_non_deleted_index];
#endif
            last_non_deleted_index--;
            total_swapped++;
          }
          if (total_swapped == actual_to_delete_count) {
            break;
          }
        }

        V[i].setOutDegree(outDegree[i] - total_swapped);
        pbbs::fetch_and_add(&numberOfSuccessfulDeletions, total_swapped);
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        inDegree[i] = V[i].getOutDegree();
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &inEdgesToDelete = currentEdgesToDelete->inEdgesToDelete;
        parallel_for(uintV j = 0;
                     j < currentEdgesToDelete->inEdgesToDelete.size(); j++) {
          uintV targetInNgh = inEdgesToDelete[j];
          uintV *currInEdges = outEdges[i];
          bool deletionSuccessful = false;
          for (intE k = 0; k < inDegree[i]; k++) {
            if (targetInNgh == currInEdges[k]) {
              bool casSuccessful;
              do {
                casSuccessful = CAS(&currInEdges[k], targetInNgh, maxValue);
              } while (currInEdges[k] != maxValue);
              if (casSuccessful) {
                deletionSuccessful = true;
                break;
              }
            }
          }
          if (deletionSuccessful == false) {
            inEdgesToDelete[j] = maxValue;
          }
        }
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        uintV *currInEdges = outEdges[i];
#ifdef EDGEDATA
        EdgeData *currInEdgeData = outEdgeData[i];
#endif

        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &inEdgesToDelete = currentEdgesToDelete->inEdgesToDelete;
        intE last_non_deleted_index;

        intE to_delete_count = inEdgesToDelete.size();

        intE actual_to_delete_count = 0;
        for (intE i = 0; i < to_delete_count; i++) {
          if (inEdgesToDelete[i] != maxValue) {
            actual_to_delete_count++;
          }
        }

        intE total_swapped = 0;

        for (intE k = inDegree[i] - 1, last_non_deleted_index = inDegree[i] - 1;
             k >= 0; k--) {
          if (currInEdges[k] == maxValue) {
            currInEdges[k] = currInEdges[last_non_deleted_index];
#ifdef EDGEDATA
            currInEdgeData[k].del();
            currInEdgeData[k] = currInEdgeData[last_non_deleted_index];
#endif
            last_non_deleted_index--;
            total_swapped++;
          }
          if (total_swapped == actual_to_delete_count) {
            break;
          }
        }
        V[i].setOutDegree(inDegree[i] - total_swapped);
        pbbs::fetch_and_add(&numberOfSuccessfulDeletions, total_swapped);
      }
    }
    intE edgeArrayIndex = 0;
    for (uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;
        vector<uintV> &inEdgesToDelete = currentEdgesToDelete->inEdgesToDelete;
#ifdef EDGEDATA
        vector<EdgeData *> &outEdgeDataToDelete =
            currentEdgesToDelete->outEdgeDataToDelete;
        vector<EdgeData *> &inEdgeDataToDelete =
            currentEdgesToDelete->inEdgeDataToDelete;
#endif

        for (intE j = 0; j < outEdgesToDelete.size(); j++) {
          if (outEdgesToDelete[j] != maxValue) {
            intE edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = outEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(outEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t" << outEdgesToDelete[j] << "\n";
            }
          }
        }

        for (intE j = 0; j < inEdgesToDelete.size(); j++) {
          if (inEdgesToDelete[j] != maxValue) {
            intE edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = inEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(inEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t" << inEdgesToDelete[j] << "\n";
            }
          }
        }
      }
    }

    free(outDegree);
    free(inDegree);
    uintE newSize = m - numberOfSuccessfulDeletions;

    m = newSize;
#ifdef EDGEDATA
    return edgeArray(ED, edgeDataWeight, edgeArrayIndex, n);
#else
    return edgeArray(ED, edgeArrayIndex, n);
#endif
  }

  edgeArray deleteEdges(edgeDeletionData &deletionsData, bool *updatedVertices,
                        bool debugFlag) {
    uintV maxValue = numeric_limits<uintV>::max();
    intE numberOfSuccessfulDeletions = 0;
    if (isSymmetric()) {
      return deleteEdges_symmetric(deletionsData, updatedVertices, debugFlag);
    }
    edge *ED = newA(edge, deletionsData.numberOfDeletions);
#ifdef EDGEDATA
    EdgeData *edgeDataWeight = newA(EdgeData, deletionsData.numberOfDeletions);
#endif
    intE *outDegree = newA(intE, n);
    intE *inDegree = newA(intE, n);

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        outDegree[i] = V[i].getOutDegree();
        inDegree[i] = V[i].getInDegree();
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;
        parallel_for(intE j = 0; j < outEdgesToDelete.size(); j++) {
          uintV targetOutNgh = outEdgesToDelete[j];
          uintV *currOutEdges = outEdges[i];

          bool deletionSuccessful = false;
          for (intE k = 0; k < outDegree[i]; k++) {
            if (targetOutNgh == currOutEdges[k]) {
              bool casSuccessful;
              do {
                casSuccessful = CAS(&currOutEdges[k], targetOutNgh, maxValue);
              } while (currOutEdges[k] != maxValue);
              if (casSuccessful) {
                deletionSuccessful = true;
                break;
              }
            }
          }
          if (deletionSuccessful == false) {
            outEdgesToDelete[j] = maxValue;
          }
        }
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &inEdgesToDelete = currentEdgesToDelete->inEdgesToDelete;
        parallel_for(intE j = 0; j < inEdgesToDelete.size(); j++) {
          uintV targetInNgh = inEdgesToDelete[j];
          uintV *currInEdges = inEdges[i];

          bool deletionSuccessful = false;
          for (intE k = 0; k < inDegree[i]; k++) {
            if (targetInNgh == currInEdges[k]) {
              bool casSuccessful;
              do {
                casSuccessful = CAS(&currInEdges[k], targetInNgh, maxValue);
              } while (currInEdges[k] != maxValue);

              if (casSuccessful) {
                deletionSuccessful = true;
                break;
              }
            }
          }
          if (deletionSuccessful == false) {
            inEdgesToDelete[j] = maxValue;
          }
        }
      }
    }

    parallel_for(uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        intE numberOfDeletions = 0;
        uintV *currOutEdges = outEdges[i];
        uintV *currInEdges = inEdges[i];
#ifdef EDGEDATA
        EdgeData *currOutEdgeData = outEdgeData[i];
        EdgeData *currInEdgeData = inEdgeData[i];
#endif

        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;
        intE last_non_deleted_index;
        intE to_delete_count = outEdgesToDelete.size();

        intE actual_to_delete_count = 0;
        for (intE i = 0; i < to_delete_count; i++) {
          if (outEdgesToDelete[i] != maxValue) {
            actual_to_delete_count++;
          }
        }
        intE total_swapped = 0;

        for (intE k = outDegree[i] - 1,
                  last_non_deleted_index = outDegree[i] - 1;
             k >= 0; k--) {
          if (currOutEdges[k] == maxValue) {
            currOutEdges[k] = currOutEdges[last_non_deleted_index];
#ifdef EDGEDATA
            currOutEdgeData[k].del();
            currOutEdgeData[k] = currOutEdgeData[last_non_deleted_index];
#endif
            last_non_deleted_index--;
            total_swapped++;
          }
          if (total_swapped == actual_to_delete_count) {
            break;
          }
        }

        V[i].setOutDegree(outDegree[i] - total_swapped);
        pbbs::fetch_and_add(&numberOfSuccessfulDeletions, total_swapped);

        vector<uintV> &inEdgesToDelete = currentEdgesToDelete->inEdgesToDelete;
        to_delete_count = inEdgesToDelete.size();
        actual_to_delete_count = 0;
        for (intE i = 0; i < to_delete_count; i++) {
          if (inEdgesToDelete[i] != maxValue) {
            actual_to_delete_count++;
          }
        }

        total_swapped = 0;

        for (intE k = inDegree[i] - 1, last_non_deleted_index = inDegree[i] - 1;
             k >= 0; k--) {
          if (currInEdges[k] == maxValue) {
            currInEdges[k] = currInEdges[last_non_deleted_index];
#ifdef EDGEDATA
            currInEdgeData[k].del();
            currInEdgeData[k] = currInEdgeData[last_non_deleted_index];
#endif
            last_non_deleted_index--;
            total_swapped++;
          }
          if (total_swapped == actual_to_delete_count) {
            break;
          }
        }
        V[i].setInDegree(inDegree[i] - total_swapped);
      }
    }
    intE edgeArrayIndex = 0;

    for (uintV i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete *currentEdgesToDelete =
            &deletionsData.getEdgeDeletionData(i);
        vector<uintV> &outEdgesToDelete =
            currentEdgesToDelete->outEdgesToDelete;
#ifdef EDGEDATA
        vector<EdgeData *> &outEdgeDataToDelete =
            currentEdgesToDelete->outEdgeDataToDelete;
#endif
        for (uintE j = 0; j < outEdgesToDelete.size(); j++) {
          if (outEdgesToDelete[j] != maxValue) {
            intE edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = outEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(outEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t" << outEdgesToDelete[j] << "\n";
            }
          }
        }
      }
    }

    free(outDegree);
    free(inDegree);
    uintE newSize = m - numberOfSuccessfulDeletions;
    m = newSize;
#ifdef EDGEDATA
    return edgeArray(ED, edgeDataWeight, numberOfSuccessfulDeletions, n);
#else
    return edgeArray(ED, edgeArrayIndex, n);
#endif
  }
  uintEE getNumberOfEdges() { return m; }

  void del() {

    if (outEdges == NULL) {
      for (uintV i = 0; i < n; i++)
        V[i].del();
    } else {
      parallel_for(uintV i = 0; i < n; i++) { 
#ifdef EDGEDATA
        parallel_for(intE j = 0; j < V[i].getOutDegree(); j++) {
          outEdgeData[i][j].del();
        }
        free(outEdgeData[i]);
#endif
        free(outEdges[i]); 
      }
      free(outEdges);
#ifdef EDGEDATA
      free(outEdgeData);
#endif
    }
    if (outEdgesArraySize != NULL) {
      free(outEdgesArraySize);
    }

    if (inEdges != NULL) {
      parallel_for(uintV i = 0; i < n; i++) { 
#ifdef EDGEDATA
        parallel_for(intE j = 0; j < V[i].getInDegree(); j++) {
          inEdgeData[i][j].del();
        }
        free(inEdgeData[i]);
#endif
        free(inEdges[i]);
      }
      free(inEdges);
#ifdef EDGEDATA
      free(inEdgeData);
#endif
    }
    free(V);

    if (inEdgesArraySize != NULL) {
      free(inEdgesArraySize);
    }

    if (outEdgeUpdates != NULL) {
      free(outEdgeUpdates);
    }
    if (inEdgeUpdates != NULL) {
      free(inEdgeUpdates);
    }
  }
};

template <class vertex> struct graph {
  vertex *V;
  uintV n;
  uintE m;
  bool transposed;
  bool symmetric;
  uintE *flags;
  Deletable *D;

  graph(vertex *_V, uintV _n, uintE _m, Deletable *_D)
      : V(_V), n(_n), m(_m), D(_D), flags(NULL), transposed(0), symmetric(0) {}

  graph(vertex *_V, uintV _n, uintE _m, Deletable *_D, uintE *_flags)
      : V(_V), n(_n), m(_m), D(_D), flags(_flags), transposed(0), symmetric(0) {
  }

  void del() {
    if (flags != NULL)
      free(flags);
    D->del();
    // free(D);
    delete D;
  }

  void setSymmetric(bool flag) {
    symmetric = flag;
    D->setSymmetric(flag);
  }

  bool isSymmetric() { return symmetric; }

  void addVertices(uintV maxVertex) {
    if (n <= maxVertex) {
      cout << "max vertex: " << maxVertex << endl;
      V = (vertex *)D->updateVertices(maxVertex);
      n = maxVertex + 1;
    }
  }

  edgeArray addEdges(edgeArray edgesToAdd, bool *updatedVertices) {
    edgeArray ea = D->addEdges(edgesToAdd, updatedVertices);
    m = D->getNumberOfEdges();
    return ea;
  }

  edgeArray deleteEdges(edgeDeletionData edgesToDelete, bool *updatedVertices,
                        bool debugFlag) {
    edgeArray ed = D->deleteEdges(edgesToDelete, updatedVertices, debugFlag);
    m = D->getNumberOfEdges();
    return ed;
  }

  void printEdges(string outputFilePath) {
#ifdef EDGEDATA
    // TODO: Add support for printing weighted edges
    return;
#endif
    cout << "M Edges: " << m << endl;
    uintV numEdgesFromDegree = 0;
    for (uintV i = 0; i < n; i++) {
      numEdgesFromDegree += V[i].getOutDegree();
    }

    if (m != numEdgesFromDegree) {
      cout << "~~~~~~~~~Edges ARE NOT EQUAL!!!!~~~~~~~~~" << endl;
      cout << "m: " << m << " NumEdges: " << numEdgesFromDegree << endl;
      abort();
    }

    ofstream outputFile;
    outputFile.open(outputFilePath, ios::out);
    outputFile << setprecision(2);
    for (uintV i = 0; i < n; ++i) {
      vertex &currentVertex = V[i];
      auto d = currentVertex.getOutDegree();
      if (d != 0) {
        insertionSort(currentVertex.getOutNeighbors(), d, ascendingF<uintV>());
      }
      for (intE j = 0; j < currentVertex.getOutDegree(); j++)
        outputFile << i << " " << currentVertex.getOutNeighbor(j) << "\n";
    }
    outputFile.close();
    outputFile.open(outputFilePath + "_inEdges", ios::out);
    outputFile << setprecision(2);
    for (uintV i = 0; i < n; ++i) {
      vertex &currentVertex = V[i];
      auto d = currentVertex.getInDegree();
      if (d != 0) {
        insertionSort(currentVertex.getInNeighbors(), d, ascendingF<uintV>());
      }
      for (intE j = 0; j < currentVertex.getInDegree(); j++) {
        outputFile << i << " " << currentVertex.getInNeighbor(j) << "\n";
      }
      // outputFile << currentVertex.getInNeighbor(j) << " " << i << "\n";
    }
    outputFile.close();
  }
};
#endif
