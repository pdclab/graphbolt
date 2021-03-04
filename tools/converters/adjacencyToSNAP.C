
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstdio>
#include <cstring>

typedef unsigned VertexType;

struct HeaderType {
  int sizeOfVertexType;
  VertexType numVertices;
  unsigned long long numEdges;
};

HeaderType header;

void readFile(std::string adjacencyListFile) {
  assert(false);
  std::ifstream alStream(adjacencyListFile);

  if(!alStream.good()) {
    std::cout << "Cannot open " << adjacencyListFile << std::endl;
    abort();
  }

  std::string line;
  while (std::getline(alStream, line)) {
    if(line[0] == '#' || line[0] == '%')
      continue;

    std::istringstream iss(line);
    VertexType dst, src;
    if (!(iss >> src >> dst)) { break; }

    if(dst == src)
      continue;

    header.numVertices = std::max(header.numVertices, std::max(src, dst)); 
    ++header.numEdges;
  }
  alStream.close();
  ++header.numVertices;
}

void readWriteFile(std::string adjacencyListFile, std::string bSFile, bool undirected, bool withheader) {
  std::ifstream alStream(adjacencyListFile);
  if(!alStream.good()) {
    std::cout << "Cannot open " << adjacencyListFile << std::endl;
    abort();
  }

  std::ofstream bSStream;
  bSStream.open(bSFile);

  if(undirected)
    header.numEdges *= 2;

  if(withheader)
    bSStream.write((char*) &header, sizeof(header));

  std::string line;
  while (std::getline(alStream, line)) {
    if(line[0] == '#' || line[0] == '%')
      continue;

    std::istringstream iss(line);
    VertexType src, count, dst;
    if (!(iss >> src >> count)) { break; }

    while(count--) {
      if(!(iss >> dst)) { break; }
      if(src == dst)
        continue;

      bSStream << src << "\t" << dst << "\n";

      if(undirected) {
        bSStream << dst << "\t" << src << "\n";
      }
    }
  }
  alStream.close();
  bSStream.close();
}

int main(int argc, char* argv[]) {
  if(argc < 4) {
    std::cout << "Please invoke like this: " << argv[0] << " --adjacencylistfile=<adjacencylistfil> --undirected=<0|1> --header=<0|1>" << std::endl;
    return -1;
  }

  std::string adjacencyListFile;
  bool undirected = false;
  bool withheader = false;

  for(int i=0; i<argc; ++i) {
    if(strncmp("--adjacencylistfile=", argv[i], 20) == 0) 
      adjacencyListFile = argv[i] + 20;
    if(strncmp("--undirected=", argv[i], 13) == 0) {
      int undir = 0;
      sscanf(argv[i] + 13, "%d", &undir);
      undirected = (undir == 0 ? false : true);
    }
    if(strncmp("--header=", argv[i], 9) == 0) {
      int hdr = 0;
      sscanf(argv[i] + 9, "%d", &hdr);
      withheader = (hdr == 0 ? false : true);
    }
  }

  if(adjacencyListFile.size() == 0) {
    std::cout << "No AdjacencyList file provided (--snapfile)." << std::endl;
    return -1; 
  } 

  std::cout << "AdjacencyList file: " << adjacencyListFile << std::endl;
  std::cout << "Unidrected: " << (undirected ? "true" : "false") << std::endl;
  std::cout << "Header: " << (withheader ? "true" : "false") << std::endl;
  std::cout << "Self-edges will be removed" << std::endl;
  std::cout << "If undirected, edge repitions might occur" << std::endl;

  header.sizeOfVertexType = sizeof(VertexType);
  header.numVertices = 0;
  header.numEdges = 0;

  if(withheader)
    readFile(adjacencyListFile);

  std::string bSFile = adjacencyListFile + ".snap";
  readWriteFile(adjacencyListFile, bSFile, undirected, withheader);

  return 0;
}
