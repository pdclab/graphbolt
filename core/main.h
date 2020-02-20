// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval
// Vora.
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
#ifndef MAIN_H
#define MAIN_H

#include "common/utils.h"
#include "graph/IO.h"
#include "graph/graph.h"
#include "graph/vertex.h"
using namespace std;

template <class vertex> void compute(graph<vertex> &, commandLine);

int parallel_main(int argc, char *argv[]) {
  commandLine P(argc, argv);
  char *iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");

  // bool useExtraBufferSpace = P.getOptionValue("-buffer");
  bool simpleFlag = P.getOptionValue("-simple");
  bool debugFlag = P.getOptionValue("-debug");

  int n_workers = P.getOptionIntValue("-nWorkers", getWorkers());
  setCustomWorkers(n_workers);

  cout << fixed;

  if (symmetric) {
    // symmetric graph
    graph<symmetricVertex> G =
        readGraph<symmetricVertex>(iFile, symmetric, simpleFlag, debugFlag);
    G.setSymmetric(true);
    cout << "Graph created" << endl;
    compute(G, P);
    G.del();
  } else {
    // asymmetric graph
    graph<asymmetricVertex> G =
        readGraph<asymmetricVertex>(iFile, symmetric, simpleFlag, debugFlag);
    cout << "Graph created" << endl;
    compute(G, P);
    G.del();
  }
}

#endif
