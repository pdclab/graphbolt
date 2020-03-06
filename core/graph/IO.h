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

#ifndef __IO_H__
#define __IO_H__

#include "graph.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;

typedef pair<uintE, uintE> intPair;
#ifdef EDGEDATA
typedef pair<uintE, pair<uintE, EdgeData *>> intWeights;
#endif

template <class E> struct pairFirstCmp {
  bool operator()(pair<uintE, E> a, pair<uintE, E> b) {
    return a.first < b.first;
  }
};

template <class E> struct getFirst {
  uintE operator()(pair<uintE, E> a) { return a.first; }
};

template <class IntType> struct pairBothCmp {
  bool operator()(pair<uintE, IntType> a, pair<uintE, IntType> b) {
    if (a.first != b.first)
      return a.first < b.first;
    return a.second < b.second;
  }
};

#ifdef EDGEDATA
struct tripleBothCmp {
  bool operator()(intWeights a, intWeights b) {
    if (a.first != b.first)
      return a.first < b.first;
    return a.second.first < b.second.first;
  }
};
#endif

struct edgeBothCmp {
  bool operator()(edge a, edge b) {
    if (a.source != b.source)
      return a.source < b.source;
    return a.destination < b.destination;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  unsigned long n; // total number of characters
  char *Chars;     // array storing all strings
  unsigned long m; // number of substrings
  char **Strings;  // pointers to strings (all should be null terminated)
  words() {}
  words(char *C, unsigned long nn, char **S, unsigned long mm)
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

_seq<char> readStringFromFile(char *fileName) {
  ifstream file(fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  unsigned long end = file.tellg();
  file.seekg(0, ios::beg);
  unsigned long n = end - file.tellg();
  char *bytes = newA(char, n + 1);
  file.read(bytes, n);
  file.close();
  return _seq<char>(bytes, n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, unsigned long n) {
  {
    parallel_for(unsigned long i = 0; i < n; i++) if (isSpace(Str[i])) Str[i] =
        0;
  }

  // mark start of words
  bool *FL = newA(bool, n);
  FL[0] = Str[0];
  {
    parallel_for(unsigned long i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1];
  }

  // offset for each start of word
  _seq<unsigned long> Off = sequence::packIndex<unsigned long>(FL, n);
  unsigned long m = Off.n;
  unsigned long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char *, m);
  { parallel_for(unsigned long j = 0; j < m; j++) SA[j] = Str + offsets[j]; }

  free(offsets);
  free(FL);
  return words(Str, n, SA, m);
}

uintE removeDuplicates(intPair *&array, uintE length, bool symmetric,
                       bool debugFlag) {
  bool *flag = newAWithZero(bool, length);
  uintE invalidCount = 0;
  parallel_for(uintE i = 1; i < length; i++) {
    if (array[i].first == array[i - 1].first &&
        array[i].second == array[i - 1].second) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].first << "\t" << array[i].second
             << "\n";
      }
    }
    if (symmetric && (array[i].first == array[i - 1].second) &&
        (array[i].second == array[i - 1].first)) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].first << "\t" << array[i].second
             << "\n";
      }
    }
  }

  uintE count = 0;
  // intPair* temp = newA(intPair, length - invalidCount);
  intPair *temp = newA(intPair, (length - invalidCount));
  for (uintE i = 0; i < length; i++) {
    if (!flag[i]) {
      temp[count].first = array[i].first;
      temp[count].second = array[i].second;
      count++;
    }
  }

  // intPair *x = array;
  free(array);
  array = temp;
  // free(array);
  free(flag);
  if (debugFlag) {
    cout << "COUNT: " << count << endl;
    cout << "INVALID COUNT: " << invalidCount << endl;
    cout << "LENGTH: " << length - invalidCount << endl;
  }
  return count;
}

#ifdef EDGEDATA
uintE removeDuplicates(intWeights *&array, uintE length, bool symmetric,
                       bool debugFlag) {
  bool *flag = newAWithZero(bool, length);
  uintE invalidCount = 0;
  parallel_for(uintE i = 1; i < length; i++) {
    if (array[i].first == array[i - 1].first &&
        array[i].second.first == array[i - 1].second.first) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].first << "\t" << array[i].second.first
             << "\n";
      }
    }
    if (symmetric && (array[i].first == array[i - 1].second.first) &&
        (array[i].second.first == array[i - 1].first)) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].first << "\t" << array[i].second.first
             << "\n";
      }
    }
  }

  uintE count = 0;
  intWeights *temp = newA(intWeights, length);
  for (uintE i = 0; i < length; i++) {
    if (!flag[i]) {
      temp[count].first = array[i].first;
      temp[count].second.first = array[i].second.first;
      temp[count].second.second = array[i].second.second;

      count++;
    }
  }

  free(array);
  array = temp;
  free(flag);
  if (debugFlag) {
    cout << "COUNT: " << count << endl;
    cout << "INVALID COUNT: " << invalidCount << endl;
    cout << "LENGTH: " << length - invalidCount << endl;
  }
  return count;
}
#endif

uintE removeDuplicates(edge *&array, uintE length, uintE maxLength,
                       bool symmetric, bool debugFlag) {
  bool *flag = newAWithZero(bool, length);
  uintE invalidCount = 0;
  parallel_for(uintE i = 1; i < length; i++) {
    if (array[i].source == array[i - 1].source &&
        array[i].destination == array[i - 1].destination) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].source << "\t" << array[i].destination
             << "\n";
      }
    }
    if (symmetric && (array[i].source == array[i - 1].destination) &&
        (array[i].destination == array[i - 1].source)) {
      flag[i] = true;
      invalidCount++;
      if (debugFlag) {
        cerr << "INVALID: " << array[i].source << "\t" << array[i].destination
             << "\n";
      }
    }
  }

  uintE count = 0;
  // edge* temp = newA(edge, length - invalidCount);
  edge *temp = newA(edge, maxLength);
  for (uintE i = 0; i < length; i++) {
    if (!flag[i]) {
      temp[count].source = array[i].source;
      temp[count].destination = array[i].destination;
#ifdef EDGEDATA
      temp[count].edgeData = array[i].edgeData;
#endif
      count++;
    }
#ifdef EDGEDATA
    else {
      array[i].edgeData->del();
    }
#endif
  }

  free(array);
  array = temp;
  free(flag);
  if (debugFlag) {
    cout << "COUNT: " << count << endl;
    cout << "INVALID COUNT: " << invalidCount << endl;
    cout << "LENGTH: " << length - invalidCount << endl;
  }
  return count;
}

template <class vertex>
graph<vertex> readGraphFromFile(char *fname, bool isSymmetric, bool simpleFlag,
                                bool debugFlag) {
  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);
#ifdef EDGEDATA
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
#else
  if (W.Strings[0] != (string) "AdjacencyGraph") {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  unsigned long len = W.m - 1;
  unsigned long n = atol(W.Strings[1]);
  unsigned long m = atol(W.Strings[2]);
  unsigned long edgeDataSize = atol(
      W.Strings[3]); // this really doesn't matter (should we even keep it?)

  // cout << "n :" << n << endl;
  // cout << "m :" << m << endl;
#ifdef EDGEDATA
  if (len != n + 2 * m + 2) {
#else
  if (len != n + m + 2) {
#endif
    cout << "Length" << len << endl;
    cout << "Bad input file" << endl;
    abort();
  }

  intE *offsets = newA(intE, n);
  uintV *edges = newA(uintV, m);
#ifdef EDGEDATA
  EdgeData *edgeData = newA(EdgeData, m);
#endif
  {
    parallel_for(unsigned long i = 0; i < n; i++) offsets[i] =
        atol(W.Strings[i + 3]);
  }
  {
    parallel_for(unsigned long i = 0; i < m; i++) {
      edges[i] = atol(W.Strings[i + n + 3]);
#ifdef EDGEDATA
      new (edgeData + i) EdgeData();
      edgeData[i].createEdgeData(W.Strings[i + n + m + 3]);
#endif
    }
  }
  W.del(); // to deal with performance bug in malloc
  vertex *v = newA(vertex, n);
  {
    parallel_for(uintV i = 0; i < n; i++) {
      intE o = offsets[i];
      intE l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      v[i].setOutDegree(l);
      v[i].setOutNeighbors(edges + o);
#ifdef EDGEDATA
      v[i].setOutEdgeDataArray(edgeData + o);
#endif
    }
  }

  intE *tOffsets;

  // TODO: ADD SYMMETRIC SUPPORT
  if (!isSymmetric) {
    tOffsets = newA(intE, n);
    { parallel_for(unsigned long i = 0; i < n; i++) tOffsets[i] = INT_E_MAX; }
#ifdef EDGEDATA
    intWeights *temp = newA(intWeights, m);
#else
    intPair *temp = newA(intPair, m);
#endif
    {
      parallel_for(unsigned long i = 0; i < n; i++) {
        intE o = offsets[i];
        for (intE j = 0; j < v[i].getOutDegree(); j++) {
#ifdef EDGEDATA
          temp[o + j] = make_pair(v[i].getOutNeighbor(j),
                                  make_pair(i, v[i].getOutEdgeData(j)));
#else
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), i);
#endif
        }
      }
    }
    // free(offsets);
#ifdef EDGEDATA
    quickSort(temp, m, tripleBothCmp());
#else
    quickSort(temp, m, pairBothCmp<intE>());
#endif
    // remove duplicates
    if (simpleFlag) {
      m = removeDuplicates(temp, m, isSymmetric, debugFlag);
    }
    tOffsets[temp[0].first] = 0;
    uintV *inEdges = newA(uintV, m);
#ifdef EDGEDATA
    inEdges[0] = temp[0].second.first;
    EdgeData *inEdgeData = newA(EdgeData, m);
    new (inEdgeData) EdgeData();
    inEdgeData[0].setEdgeDataFromPtr(temp[0].second.second);
#else
    inEdges[0] = temp[0].second;
#endif
    {
      parallel_for(unsigned long i = 1; i < m; i++) {
#ifdef EDGEDATA
        inEdges[i] = temp[i].second.first;
        new (inEdgeData + i) EdgeData();
        inEdgeData[i].setEdgeDataFromPtr(temp[i].second.second);
#else
        inEdges[i] = temp[i].second;
#endif
        if (temp[i].first != temp[i - 1].first) {
          tOffsets[temp[i].first] = i;
        }
      }
    }

    free(temp);

    // fill in offsets of degree 0 vertices by taking closest non-zero
    // offset to the right
    sequence::scanIBack(tOffsets, tOffsets, n, minF<intE>(), (intE)m);

    {
      parallel_for(unsigned long i = 0; i < n; i++) {
        uintE o = tOffsets[i];
        uintE l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        v[i].setInDegree(l);
        v[i].setInNeighbors(inEdges + o);
#ifdef EDGEDATA
        v[i].setInEdgeDataArray(inEdgeData + o);
#endif
      }
    }

#ifdef EDGEDATA
    AdjacencyRep<vertex> *mem = new AdjacencyRep<vertex>(
        v, n, m, edges, inEdges, offsets, tOffsets, edgeData, inEdgeData);
#else
    AdjacencyRep<vertex> *mem =
        new AdjacencyRep<vertex>(v, n, m, edges, inEdges, offsets, tOffsets);
#endif
    return graph<vertex>(v, n, m, mem);
  } else {
#ifdef EDGEDATA
    AdjacencyRep<vertex> *mem = new AdjacencyRep<vertex>(
        v, n, m, edges, NULL, offsets, NULL, edgeData, NULL);
#else
    AdjacencyRep<vertex> *mem =
        new AdjacencyRep<vertex>(v, n, m, edges, NULL, offsets, NULL);
#endif
    return graph<vertex>(v, n, m, mem);
  }
}

#ifdef EDGEDATA
template <class vertex>
void printWeightedGraph(string outputFilePath, graph<vertex> G) {
  cout << "Weighted graph" << outputFilePath << endl;
  intWeights *pairedEdges = newA(intWeights, G.m);
  cout << "M Edges: " << G.m << endl;
  intE numEdgesFromDegree = 0;
  for (uintV i = 0; i < G.n; i++) {
    numEdgesFromDegree += G.V[i].getOutDegree();
  }

  if (G.m != numEdgesFromDegree) {
    cout << "~~~~~~~~~Edges ARE NOT EQUAL!!!!~~~~~~~~~" << endl;
    cout << "G.m: " << G.m << " NumEdges: " << numEdgesFromDegree << endl;
    abort();
  }

  intE offset = 0;
  for (uintV i = 0; i < G.n; i++) {
    for (intE j = 0; j < G.V[i].getOutDegree(); j++) {
      // file << i << " " << G.V[i].Neighbors[j] << " " << G.V[i].nghWeights[j]
      // << "\n"; EdgeData *w = G.V[i].getOutEdgeData(j)->copyEdgeData();
      pairedEdges[offset + j] = make_pair(
          i, make_pair(G.V[i].getOutNeighbor(j), G.V[i].getOutEdgeData(j)));
    }
    offset += G.V[i].getOutDegree();
  }
  cout << "Offset: " << offset << endl;
  quickSort(pairedEdges, offset, tripleBothCmp());

  ofstream outputFile;
  outputFile.open(outputFilePath, ios::out);
  outputFile << setprecision(2);

  for (uintE i = 0; i < offset; i++) {
    std::string weight = (pairedEdges[i].second.second)->print();
    outputFile << pairedEdges[i].first << " " << pairedEdges[i].second.first
               << " " << weight << endl;
    // outputFile << pairedEdges[i].first << " " << pairedEdges[i].second.first
    // << " " << pairedEdges[i].second.second->print() << endl;
  }
  outputFile.close();
  free(pairedEdges);
}
#endif

template <class vertex>
graph<vertex> readGraph(char *iFile, bool symmetric, bool isSimple,
                        bool debugFlag) {
  return readGraphFromFile<vertex>(iFile, symmetric, isSimple, debugFlag);
}
#endif
