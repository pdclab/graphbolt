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
  long source;
  long destination;
  edge(long f, long s) : source(f), destination(s) {}
};

struct edgeArray {
  edge *E;
  unsigned long size;
  unsigned long maxVertex;
  void del() {
    if (E != nullptr)
      free(E);
  }
  edgeArray() : size(0) { E = nullptr; }
  edgeArray(edge *EE, unsigned long _size, unsigned long _maxVertex)
      : E(EE), size(_size), maxVertex(_maxVertex) {}
};

struct edgesToDelete {
  vector<uintT> outEdgesToDelete;
  vector<uintT> inEdgesToDelete;

  edgesToDelete() {}
  edgesToDelete(vector<uintT> _outEdgesToDelete, vector<uintT> _inEdgesToDelete)
      : outEdgesToDelete(_outEdgesToDelete), inEdgesToDelete(_inEdgesToDelete) {
  }

  inline bool hasOutNeighbor(uintT outNeighbor) {
    vector<uintT>::iterator index =
        find(outEdgesToDelete.begin(), outEdgesToDelete.end(), outNeighbor);
    if (index != outEdgesToDelete.end()) {
      return true;
    }
    return false;
  }

  inline bool hasInNeighbor(uintT inNeighbor) {
    vector<uintT>::iterator index =
        find(inEdgesToDelete.begin(), inEdgesToDelete.end(), inNeighbor);
    if (index != inEdgesToDelete.end()) {
      return true;
    }
    return false;
  }

  inline long getIndexOfOutNeighbor(uintT outNeighbor, long start = 0) {

    vector<uintT>::iterator it = find(outEdgesToDelete.begin() + start,
                                      outEdgesToDelete.end(), outNeighbor);
    if (it != outEdgesToDelete.end()) {
      return it - outEdgesToDelete.begin();
    }
    return -1;
  }

  inline long getIndexOfInNeighbor(uintT inNeighbor, long start = 0) {
    vector<uintT>::iterator it = find(inEdgesToDelete.begin() + start,
                                      inEdgesToDelete.end(), inNeighbor);
    if (it != inEdgesToDelete.end()) {
      return it - inEdgesToDelete.begin();
    }
    return -1;
  }

  inline void clear() {
    outEdgesToDelete.clear();
    inEdgesToDelete.clear();
  }

  inline void insertOutEdge(uintT v) { outEdgesToDelete.push_back(v); }

  inline void insertInEdge(uintT v) { inEdgesToDelete.push_back(v); }
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
  uintT n;
  bool *updatedVertices;

  edgeDeletionData(uintT _n) : n(_n) {
    updatedVertices = newA(bool, n);
    dataMap = new edgesToDelete[n];
    parallel_for(uintT i = 0; i < n; i++) { updatedVertices[i] = 0; }
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

    parallel_for(intT i = 0; i < numberOfDeletions; i++) {
      intT prev = (i == 0) ? 0 : i - 1;
      if ((i == 0) || (edgesArray[prev].source != edgesArray[i].source)) {
        intT j = i;
        updatedVertices[edgesArray[i].source] = 1;
        while ((j < numberOfDeletions) &&
               (edgesArray[i].source == edgesArray[j].source)) {
          dataMap[edgesArray[j].source].insertOutEdge(
              edgesArray[j].destination);
          j++;
        }
      }
    }

    // sort on dest
    quickSort(edgesArray, numberOfDeletions, MinDesCmp<edge>());
    parallel_for(intT i = 0; i < numberOfDeletions; i++) {
      intT prev = (i == 0) ? 0 : i - 1;
      if (i == 0 || edgesArray[prev].destination != edgesArray[i].destination) {
        intT j = i;
        updatedVertices[edgesArray[i].destination] = 1;
        while (j < numberOfDeletions &&
               edgesArray[i].destination == edgesArray[j].destination) {
          dataMap[edgesArray[j].destination].insertInEdge(edgesArray[j].source);
          j++;
        }
      }
    }
  }

  edgesToDelete &getEdgeDeletionData(uintT vertex) { return dataMap[vertex]; }

  void reset() {
    numberOfDeletions = 0;
    parallel_for(uintT i = 0; i < n; i++) {
      if (updatedVertices[i] == 1) {
        // Delete the entries in dataMap corresponding to the vertex i
        dataMap[i].clear();
        updatedVertices[i] = 0;
      }
    }
  }

  void updateNumVertices(intT maxVertex) {
    if (n) {
      delete[] dataMap;
    }
    updatedVertices = renewA(bool, updatedVertices, maxVertex);
    parallel_for(uintT i = 0; i < maxVertex; i++) { updatedVertices[i] = 0; }
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
  virtual void *updateVertices(uintT _verticesSize) = 0;
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
  uintE **outEdges = NULL, **inEdges = NULL;
  uintE *outEdgesArraySize = NULL, *inEdgesArraySize = NULL;
  // Indicates whether the graph is directed or undirected. TRUE -> indicates
  // undirected graph
  bool symmetric;

  uintE *outEdgeArrayLengths;
  uintE *inEdgeArrayLengths;

  uintT *inEdgeUpdates = NULL;
  uintT *outEdgeUpdates = NULL;

  AdjacencyRep(vertex *VV, unsigned long nn, unsigned long mm, uintE *ai,
                      uintE *_inEdges = NULL, uintE *_outEdgeOffsets = NULL,
                      uintE *_inEdgeOffsets = NULL)
      : V(VV), n(nn), m(mm), symmetric(false) {
    outEdges = newA(uintE *, nn);
    outEdgesArraySize = newA(uintE, nn);
    outEdgeUpdates = newA(uintT, nn);

#ifdef EDGEDATA
    if (_outEdgeData != NULL) {
      outEdgeData = newA(EdgeData *, nn);
    }
#endif

    if (_inEdges != NULL) {
      inEdges = newA(uintE *, nn);
      inEdgesArraySize = newA(uintE, nn);
      inEdgeUpdates = newA(uintT, nn);
#ifdef EDGEDATA
      inEdgeData = newA(EdgeData *, nn);
#endif
    }

    parallel_for(unsigned long i = 0; i < nn; i++) {
      uintE outEdgesSize = V[i].getOutDegree();
      outEdges[i] = newA(uintE, outEdgesSize);
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
        inEdges[i] = newA(uintE, inEdgesSize);
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
      parallel_for(uintT i = 0; i < mm; i++) { _outEdgeData[i].del(); }
      free(_outEdgeData);
#endif
    }

    if (_inEdges != NULL) {
      free(_inEdges);
      free(_inEdgeOffsets);
#ifdef EDGEDATA
      parallel_for(uintT i = 0; i < mm; i++) { _inEdgeData[i].del(); }
      free(_inEdgeData);
#endif
    }
  }

  // template <class intT>
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

  void *updateVertices(uintT maxVertex) {
    if (isSymmetric()) {
      return updateVertices_symmetric(maxVertex);
    }

    if (maxVertex >= n) {
      uintT currentVertexSize = n;
      n = maxVertex + 1;
      V = renewA(vertex, V, n);
      outEdges = renewA(uintE *, outEdges, n);
      inEdges = renewA(uintE *, inEdges, n);
      outEdgesArraySize = renewA(uintE, outEdgesArraySize, n);
      inEdgesArraySize = renewA(uintE, inEdgesArraySize, n);
      outEdgeUpdates = renewA(uintT, outEdgeUpdates, n);
      inEdgeUpdates = renewA(uintT, inEdgeUpdates, n);
#ifdef EDGEDATA
      outEdgeData = renewA(EdgeData *, outEdgeData, n);
      inEdgeData = renewA(EdgeData *, inEdgeData, n);
#endif

      parallel_for(uintT i = 0; i < currentVertexSize; i++) {
        V[i].setOutNeighbors(outEdges[i]);
        V[i].setInNeighbors(inEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
        V[i].setInEdgeDataArray(inEdgeData[i]);
#endif
      }

      parallel_for(uintT i = currentVertexSize; i < n; i++) {
        V[i].setOutDegree(0);
        V[i].setInDegree(0);
        outEdges[i] = newA(uintE, 0);
        inEdges[i] = newA(uintE, 0);
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

  void *updateVertices_symmetric(uintT maxVertex) {
    if (maxVertex > n) {
      uintT currentVertexSize = n;
      n = maxVertex + 1;
      V = renewA(vertex, V, n);
      outEdges = renewA(uintE *, outEdges, n);
      outEdgesArraySize = renewA(uintE, outEdgesArraySize, n);
      outEdgeUpdates = renewA(uintT, outEdgeUpdates, n);
#ifdef EDGEDATA
      outEdgeData = renewA(EdgeData *, outEdgeData, n);
#endif

      parallel_for(uintT i = 0; i < currentVertexSize; i++) {
        V[i].setOutNeighbors(outEdges[i]);
#ifdef EDGEDATA
        V[i].setOutEdgeDataArray(outEdgeData[i]);
#endif
      }

      parallel_for(uintT i = currentVertexSize; i < n; i++) {
        V[i].setOutDegree(0);
        outEdges[i] = newA(uintE, 0);
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
    parallel_for(uintE i = 0; i < n; i++) { outEdgeUpdates[i] = 0; }

    edge *E = edgesToAdd.E;
    uintE size = edgesToAdd.size;
    // for each new edge, we increment the outDegreeOffset of the source
    // and we increment the inDegreeOffset of the destination
    parallel_for(uintE i = 0; i < size; i++) {
      uint source = E[i].source;
      uint destination = E[i].destination;

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
        outEdges[i] = renewA(uintE, outEdges[i], outEdgesArraySize[i]);
        V[i].setOutNeighbors(outEdges[i]);
      }
#else
      if (outEdgeUpdates[i]) {
        outEdges[i] = renewA(uintE, outEdges[i],
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
      uintE source = E[i].source;
      uintE destination = E[i].destination;
#ifdef EDGEDATA
      EdgeData *edgeData = E[i].edgeData;
#endif
      vertex &sourceVertex = V[source];
      vertex &destinationVertex = V[destination];
      int outIndex = sourceVertex.fetchAndAddOutDegree(1);
      outEdges[source][outIndex] = destination;
      int inIndex = destinationVertex.fetchAndAddOutDegree(1);
      outEdges[destination][inIndex] = source;
#ifdef EDGEDATA
      new (outEdgeData[source] + outIndex) EdgeData();
      outEdgeData[source][outIndex].setEdgeDataFromPtr(edgeData);
      new (outEdgeData[destination] + inIndex) EdgeData();
      outEdgeData[destination][inIndex].setEdgeDataFromPtr(edgeData);
#endif
    }

    long curr_size = edgesToAdd.size;
    edge *symE = newA(edge, curr_size * 2);
    long j = 0;
    for (long i = 0; i < curr_size; i++) {
      symE[j].source = edgesToAdd.E[i].source;
      symE[j].destination = edgesToAdd.E[i].destination;
      j++;
      symE[j].source = edgesToAdd.E[i].destination;
      symE[j].destination = edgesToAdd.E[i].source;
      j++;
    }
    long max_vertex = edgesToAdd.maxVertex;
    edgesToAdd.del();
    edgesToAdd = edgeArray(symE, j, max_vertex);

    uintE newSize = m + (edgesToAdd.size);
    m = newSize;
    return edgesToAdd;
  }

  edgeArray addEdges(edgeArray &edgesToAdd, bool *updatedVertices) {
    if (isSymmetric()) {
      return addEdges_symmetric(edgesToAdd, updatedVertices);
    }

    parallel_for(uintT i = 0; i < n; i++) {
      outEdgeUpdates[i] = 0;
      inEdgeUpdates[i] = 0;
    }

    edge *E = edgesToAdd.E;
    uint size = edgesToAdd.size;
    // for each new edge, we increment the outDegreeOffset of the source
    // and we increment the inDegreeOffset of the destination
    parallel_for(uintE i = 0; i < size; i++) {
      uintT source = E[i].source;
      uintT destination = E[i].destination;

      pbbs::write_add(&outEdgeUpdates[source], 1);

      pbbs::write_add(&inEdgeUpdates[destination], 1);

      updatedVertices[source] = 1;
      updatedVertices[destination] = 1;
    }

    parallel_for(uintE i = 0; i < n; i++) {
#ifdef INCLUDEEXTRABUFFERSPACE
      if (outEdgeUpdates[i] &&
          (outEdgeUpdates[i] + V[i].getOutDegree()) > outEdgesArraySize[i]) {
        outEdgesArraySize[i] =
            V[i].getOutDegree() + outEdgeUpdates[i] + EDGE_ARRAY_BUFFER;
        outEdges[i] = renewA(uintE, outEdges[i], outEdgesArraySize[i]);
        V[i].setOutNeighbors(outEdges[i]);
      }
      if (inEdgeUpdates[i] &&
          (inEdgeUpdates[i] + V[i].getInDegree()) > inEdgesArraySize[i]) {
        inEdgesArraySize[i] =
            V[i].getInDegree() + inEdgeUpdates[i] + EDGE_ARRAY_BUFFER;
        inEdges[i] = renewA(uintE, inEdges[i], inEdgesArraySize[i]);
        V[i].setInNeighbors(inEdges[i]);
      }
#else
      if (outEdgeUpdates[i]) {
        outEdges[i] = renewA(uintE, outEdges[i],
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
            renewA(uintE, inEdges[i], (V[i].getInDegree() + inEdgeUpdates[i]));
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
      uintE source = E[i].source;
      uintE destination = E[i].destination;
#ifdef EDGEDATA
      EdgeData *edgeData = E[i].edgeData;
#endif
      vertex &sourceVertex = V[source];
      vertex &destinationVertex = V[destination];

      int outIndex = sourceVertex.fetchAndAddOutDegree(1);
      outEdges[source][outIndex] = destination;
      int inIndex = destinationVertex.fetchAndAddInDegree(1);
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
    uintE maxValue = numeric_limits<uintE>::max();
    intT numberOfSuccessfulDeletions = 0;

    edge *ED = newA(edge, deletionsData.numberOfDeletions * 2);
#ifdef EDGEDATA
    EdgeData *edgeDataWeight = newA(EdgeData, deletionsData.numberOfDeletions);
#endif
    intT edgeArrayIndex = 0;

    intT **outIndicesArray = newA(intT *, n);
    intT **inIndicesArray = newA(intT *, n);

    intT *outIndex = newA(intT, n);
    intT *inIndex = newA(intT, n);

    intT *outDegree = newA(intT, n);
    intT *inDegree = newA(intT, n);

    bool **outFlagArray = newA(bool *, n);
    bool **inFlagArray = newA(bool *, n);

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        outDegree[i] = V[i].getOutDegree();

        if (currentEdgesToDelete.outEdgesToDelete.size() > 0) {
          outIndicesArray[i] =
              newA(intT, currentEdgesToDelete.outEdgesToDelete.size());
          outFlagArray[i] =
              newA(bool, currentEdgesToDelete.outEdgesToDelete.size());
        }
        if (currentEdgesToDelete.inEdgesToDelete.size() > 0) {
          inIndicesArray[i] =
              newA(intT, currentEdgesToDelete.inEdgesToDelete.size());
          inFlagArray[i] =
              newA(bool, currentEdgesToDelete.inEdgesToDelete.size());
        }
        outIndex[i] = 0;
        inIndex[i] = 0;
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *outFlag = outFlagArray[i];
        parallel_for(intT j = 0;
                     j < currentEdgesToDelete.outEdgesToDelete.size(); j++) {
          outFlag[j] = false;
          uintT targetOutNgh = currentEdgesToDelete.outEdgesToDelete[j];
          uintE *currOutEdges = outEdges[i];
  
          for (intT k = 0; k < outDegree[i]; k++) {
            if (targetOutNgh == currOutEdges[k]) {
              bool casSuccessful =
                  CAS(&currOutEdges[k], targetOutNgh, maxValue);
              if (casSuccessful) {
                outIndicesArray[i][j] = k;
                outFlag[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        intT *outIndices = outIndicesArray[i];
        bool *currOutFlags = outFlagArray[i];
        uintE *currOutEdges = outEdges[i];

        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        intE last_non_deleted_index;
        long to_delete_count = currentEdgesToDelete.outEdgesToDelete.size();
        long actual_to_delete_count = 0;
        for (long i = 0; i < to_delete_count; i++) {
          if (currOutFlags[i] == true) {
            actual_to_delete_count++;
          }
        }
        int total_swapped = 0;

        for (intE k = outDegree[i] - 1,
                  last_non_deleted_index = outDegree[i] - 1;
             k >= 0; k--) {
          if (currOutEdges[k] == maxValue) {
            currOutEdges[k] = currOutEdges[last_non_deleted_index];
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

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        inDegree[i] = V[i].getOutDegree();
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *inFlag = inFlagArray[i];
        parallel_for(intT j = 0;
                     j < currentEdgesToDelete.inEdgesToDelete.size(); j++) {
          inFlag[j] = false;
          uintT targetInNgh = currentEdgesToDelete.inEdgesToDelete[j];
          uintE *currInEdges = outEdges[i];
          intT *currInIndexAddr = &inIndex[i];

          for (intT k = 0; k < inDegree[i]; k++) {
            if (targetInNgh == currInEdges[k]) {
              bool casSuccessful = CAS(&currInEdges[k], targetInNgh, maxValue);
              if (casSuccessful) {
                inIndicesArray[i][j] = k;
                inFlag[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        intT *inIndices = inIndicesArray[i];
        bool *currInFlags = inFlagArray[i];
        uintE *currInEdges = outEdges[i];

        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        intE last_non_deleted_index;

        long to_delete_count = currentEdgesToDelete.inEdgesToDelete.size();
        long actual_to_delete_count = 0;
        for (long i = 0; i < to_delete_count; i++) {
          if (currInFlags[i] == 1) {
            actual_to_delete_count++;
          }
        }

        int total_swapped = 0;

        // TODO: Add EdgeData Back
        for (intE k = inDegree[i] - 1, last_non_deleted_index = inDegree[i] - 1;
             k >= 0; k--) {
          if (currInEdges[k] == maxValue) {
            currInEdges[k] = currInEdges[last_non_deleted_index];
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

    for (intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *outFlag = outFlagArray[i];
        bool *inFlag = inFlagArray[i];
        for (uintE j = 0; j < currentEdgesToDelete.outEdgesToDelete.size();
             j++) {
          if (outFlag[j] == true) {
            intT edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = currentEdgesToDelete.outEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(
                currentEdgesToDelete.outEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t"
                   << currentEdgesToDelete.outEdgesToDelete[j] << "\n";
            }
          }
        }

        for (uintE j = 0; j < currentEdgesToDelete.inEdgesToDelete.size();
             j++) {
          if (inFlag[j] == true) {
            intT edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = currentEdgesToDelete.inEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(
                currentEdgesToDelete.outEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t"
                   << currentEdgesToDelete.outEdgesToDelete[j] << "\n";
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        if (deletionsData.getEdgeDeletionData(i).outEdgesToDelete.size() > 0) {
          free(outIndicesArray[i]);
          free(outFlagArray[i]);
        }
        if (deletionsData.getEdgeDeletionData(i).inEdgesToDelete.size() > 0) {
          free(inIndicesArray[i]);
          free(inFlagArray[i]);
        }
      }
    }
    free(outIndicesArray);
    free(inIndicesArray);
    free(outFlagArray);
    free(inFlagArray);
    free(outIndex);
    free(inIndex);
    free(outDegree);
    free(inDegree);
    uintE newSize = m - numberOfSuccessfulDeletions;

    m = newSize;
#ifndef EDGEDATA
    return edgeArray(ED, edgeArrayIndex, n);
#else
    return edgeArray(ED, edgeDataWeight, edgeArrayIndex, n);
#endif
  }

  edgeArray deleteEdges(edgeDeletionData &deletionsData, bool *updatedVertices,
                        bool debugFlag) {
    uintE maxValue = numeric_limits<uintE>::max();
    intT numberOfSuccessfulDeletions = 0;
    if (isSymmetric()) {
      return deleteEdges_symmetric(deletionsData, updatedVertices, debugFlag);
    }
    edge *ED = newA(edge, deletionsData.numberOfDeletions);
#ifdef EDGEDATA
    EdgeData *edgeDataWeight = newA(EdgeData, deletionsData.numberOfDeletions);
#endif
    intT edgeArrayIndex = 0;

    intT **outIndicesArray = newA(intT *, n);
    intT **inIndicesArray = newA(intT *, n);

    intT *outIndex = newA(intT, n);
    intT *inIndex = newA(intT, n);

    intT *outDegree = newA(intT, n);
    intT *inDegree = newA(intT, n);

    bool **outFlagArray = newA(bool *, n);
    bool **inFlagArray = newA(bool *, n);

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        outDegree[i] = V[i].getOutDegree();
        inDegree[i] = V[i].getInDegree();

        if (currentEdgesToDelete.outEdgesToDelete.size() > 0) {
          outIndicesArray[i] =
              newA(intT, currentEdgesToDelete.outEdgesToDelete.size());
          outFlagArray[i] =
              newA(bool, currentEdgesToDelete.outEdgesToDelete.size());
        }
        if (currentEdgesToDelete.inEdgesToDelete.size() > 0) {
          inIndicesArray[i] =
              newA(intT, currentEdgesToDelete.inEdgesToDelete.size());
          inFlagArray[i] =
              newA(bool, currentEdgesToDelete.inEdgesToDelete.size());
        }
        outIndex[i] = 0;
        inIndex[i] = 0;
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *outFlag = outFlagArray[i];
        parallel_for(intT j = 0;
                     j < currentEdgesToDelete.outEdgesToDelete.size(); j++) {
          outFlag[j] = false;
          uintT targetOutNgh = currentEdgesToDelete.outEdgesToDelete[j];
          uintE *currOutEdges = outEdges[i];
  
          for (intT k = 0; k < outDegree[i]; k++) {
            if (targetOutNgh == currOutEdges[k]) {
              bool casSuccessful =
                  CAS(&currOutEdges[k], targetOutNgh, maxValue);
              if (casSuccessful) {
                outIndicesArray[i][j] = k;
                outFlag[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *inFlag = inFlagArray[i];
        parallel_for(intT j = 0;
                     j < currentEdgesToDelete.inEdgesToDelete.size(); j++) {
          inFlag[j] = false;
          uintT targetInNgh = currentEdgesToDelete.inEdgesToDelete[j];
          uintE *currInEdges = inEdges[i];
          intT *currInIndexAddr = &inIndex[i];

          for (intT k = 0; k < inDegree[i]; k++) {
            if (targetInNgh == currInEdges[k]) {
              bool casSuccessful = CAS(&currInEdges[k], targetInNgh, maxValue);
              if (casSuccessful) {
                inIndicesArray[i][j] = k;
                inFlag[j] = true;
                break;
              }
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        intT numberOfDeletions = 0;
        intT *outIndices = outIndicesArray[i];
        intT *inIndices = inIndicesArray[i];
        bool *currOutFlags = outFlagArray[i];
        bool *currInFlags = inFlagArray[i];
        uintE *currOutEdges = outEdges[i];
        uintE *currInEdges = inEdges[i];

        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        intE last_non_deleted_index;
        long to_delete_count = currentEdgesToDelete.outEdgesToDelete.size();
        long actual_to_delete_count = 0;
        for (long i = 0; i < to_delete_count; i++) {
          if (currOutFlags[i] == true) {
            actual_to_delete_count++;
          }
        }
        int total_swapped = 0;

        for (intE k = outDegree[i] - 1,
                  last_non_deleted_index = outDegree[i] - 1;
             k >= 0; k--) {
          if (currOutEdges[k] == maxValue) {
            currOutEdges[k] = currOutEdges[last_non_deleted_index];
            last_non_deleted_index--;
            total_swapped++;
          }
          if (total_swapped == actual_to_delete_count) {
            break;
          }
        }

        V[i].setOutDegree(outDegree[i] - total_swapped);
        pbbs::fetch_and_add(&numberOfSuccessfulDeletions, total_swapped);

        to_delete_count = currentEdgesToDelete.inEdgesToDelete.size();
        actual_to_delete_count = 0;
        for (long i = 0; i < to_delete_count; i++) {
          if (currInFlags[i] == 1) {
            actual_to_delete_count++;
          }
        }

        total_swapped = 0;

        // TODO: Add EdgeData Back
        for (intE k = inDegree[i] - 1, last_non_deleted_index = inDegree[i] - 1;
             k >= 0; k--) {
          if (currInEdges[k] == maxValue) {
            currInEdges[k] = currInEdges[last_non_deleted_index];
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

    for (intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        edgesToDelete currentEdgesToDelete =
            deletionsData.getEdgeDeletionData(i);
        bool *outFlag = outFlagArray[i];
        for (uintE j = 0; j < currentEdgesToDelete.outEdgesToDelete.size();
             j++) {
          if (outFlag[j] == true) {
            intT edIndex = edgeArrayIndex++;
            ED[edIndex].source = i;
            ED[edIndex].destination = currentEdgesToDelete.outEdgesToDelete[j];
#ifdef EDGEDATA
            new (edgeDataWeight + edIndex) EdgeData();
            ED[edIndex].edgeData = &edgeDataWeight[edIndex];
            ED[edIndex].edgeData->setEdgeDataFromPtr(
                currentEdgesToDelete.outEdgeDataToDelete[j]);
#endif
          } else {
            if (debugFlag) {
              cerr << "INVALID: " << i << "\t"
                   << currentEdgesToDelete.outEdgesToDelete[j] << "\n";
            }
          }
        }
      }
    }

    parallel_for(intT i = 0; i < n; i++) {
      if (deletionsData.updatedVertices[i] == 1) {
        if (deletionsData.getEdgeDeletionData(i).outEdgesToDelete.size() > 0) {
          free(outIndicesArray[i]);
          free(outFlagArray[i]);
        }
        if (deletionsData.getEdgeDeletionData(i).inEdgesToDelete.size() > 0) {
          free(inIndicesArray[i]);
          free(inFlagArray[i]);
        }
      }
    }
    free(outIndicesArray);
    free(inIndicesArray);
    free(outFlagArray);
    free(inFlagArray);
    free(outIndex);
    free(inIndex);
    free(outDegree);
    free(inDegree);
    uintE newSize = m - numberOfSuccessfulDeletions;
    m = newSize;
#ifndef EDGEDATA
    return edgeArray(ED, edgeArrayIndex, n);
#else
    return edgeArray(ED, edgeDataWeight, numberOfSuccessfulDeletions, n);
#endif
  }
  uintEE getNumberOfEdges() { return m; }

  void del() {

    if (outEdges == NULL) {
      for (unsigned long i = 0; i < n; i++)
        V[i].del();
    } else {
      parallel_for(unsigned long i = 0; i < n; i++) { free(outEdges[i]); }
      free(outEdges);
    }
    if (outEdgesArraySize != NULL) {
      free(outEdgesArraySize);
    }

    free(V);
    if (inEdges != NULL) {
      parallel_for(unsigned long i = 0; i < n; i++) { free(inEdges[i]); }
      free(inEdges);
    }

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
  unsigned long n;
  unsigned long m;
  bool transposed;
  bool symmetric;
  uintE *flags;
  Deletable *D;

  graph(vertex *_V, unsigned long _n, unsigned long _m, Deletable *_D)
      : V(_V), n(_n), m(_m), D(_D), flags(NULL), transposed(0), symmetric(0) {}

  graph(vertex *_V, unsigned long _n, unsigned long _m, Deletable *_D,
        uintE *_flags)
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

  void addVertices(uintT maxVertex) {
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
    cout << "M Edges: " << m << endl;
    uintT numEdgesFromDegree = 0;
    for (uintT i = 0; i < n; i++) {
      numEdgesFromDegree += V[i].getOutDegree();
    }

    if (m != numEdgesFromDegree) {
      cout << "~~~~~~~~~Edges ARE NOT EQUAL!!!!~~~~~~~~~" << endl;
      cout << "m: " << m << " NumEdges: " << numEdgesFromDegree << endl;
      abort();
    }

    parallel_for(intT i = 0; i < n; i++) {
      quickSort(V[i].getOutNeighbors(), V[i].getOutDegree(), SimpleCmp<intT>());
    }

    parallel_for(intT i = 0; i < n; i++) {
      quickSort(V[i].getInNeighbors(), V[i].getInDegree(), SimpleCmp<intT>());
    }

    ofstream outputFile;
    outputFile.open(outputFilePath, ios::out);
    outputFile << setprecision(2);
    for (uintT i = 0; i < n; ++i) {
      vertex &currentVertex = V[i];
      auto d = currentVertex.getOutDegree();
      if (d != 0) {
        insertionSort(currentVertex.getOutNeighbors(), d, ascendingF<uintT>());
      }
      for (uintE j = 0; j < currentVertex.getOutDegree(); j++)
        outputFile << i << " " << currentVertex.getOutNeighbor(j) << "\n";
    }
    outputFile.close();
    outputFile.open(outputFilePath + "_inEdges", ios::out);
    outputFile << setprecision(2);
    for (uintT i = 0; i < n; ++i) {
      vertex &currentVertex = V[i];
      auto d = currentVertex.getInDegree();
      if (d != 0) {
        insertionSort(currentVertex.getInNeighbors(), d, ascendingF<uintT>());
      }
      for (uintE j = 0; j < currentVertex.getInDegree(); j++) {
        outputFile << i << " " << currentVertex.getInNeighbor(j) << "\n";
      }
      // outputFile << currentVertex.getInNeighbor(j) << " " << i << "\n";
    }
    outputFile.close();
  }
};
#endif
