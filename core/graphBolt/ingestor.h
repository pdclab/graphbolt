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

#ifndef INGESTOR_H
#define INGESTOR_H

#include "../common/parseCommandLine.h"
#include "../common/utils.h"
#include "../graph/IO.h"
#include "../graph/graph.h"
#include <string>

/**
 * Used to extract values from the binary stream
 **/
struct StreamEdge {
  uintV source;
  uintV destination;
  char edgeType;
};

// TODO : Determine what should be abstracted and stuff
template <class vertex> class Ingestor {

public:
  commandLine config;
  graph<vertex> &my_graph;
  long n;
  bool is_symmetric;
  edgeArray edge_additions;
  edgeArray edge_deletions;
  edgeArray edge_deletions_temp;
  edgeDeletionData deletions_data;
  bool *updated_vertices;

  ifstream stream_file;

  long max_batch_size;

  long number_of_batches;
  long current_batch;
  char *stream_path;

  bool simple_flag;
  bool fixed_batch_flag;
  bool enforce_edge_validity_flag;
  bool debug_flag;
  bool stream_closed = false;

  Ingestor(graph<vertex> &_my_graph, commandLine _config)
      : my_graph(_my_graph), config(_config), n(_my_graph.n),
        deletions_data(_my_graph.n), current_batch(0) {

    updated_vertices = newA(bool, n);

    is_symmetric = config.getOptionValue("-s");
    stream_path = config.getOptionValue("-streamPath");
    simple_flag = config.getOptionValue("-simple");
    fixed_batch_flag = config.getOptionValue("-fixedBatchSize");
    enforce_edge_validity_flag = config.getOptionValue("-enforceEdgeValidity");
    debug_flag = config.getOptionValue("-debug");
    max_batch_size = config.getOptionLongValue("-nEdges", 0);
    if (max_batch_size == 0) {
      std::cout
          << "WARNING: Max Batch Size = 0, Please assign a value to '-nEdges'"
          << std::endl;
    }

    number_of_batches = config.getOptionLongValue("-numberOfUpdateBatches", 1);
    cout << "Number of batches: " << number_of_batches << endl;
  }

  ~Ingestor() {
    if (current_batch > number_of_batches) {
      edge_deletions.del();
    } else {
      edge_deletions_temp.del();
    }
    edge_additions.del();
    stream_file.close();
    deletions_data.del();
    if (n > 0)
      free(updated_vertices);
  }

  edgeArray &getEdgeAdditions() {
    cout << "Edge Additions in batch: " << edge_additions.size << endl;
    return edge_additions;
  }

  edgeArray &getEdgeDeletions() {
    cout << "Edge Deletions in batch: " << edge_deletions.size << endl;
    return edge_deletions;
  }

  void cleanup() {
    cout << "Current_batch: " << current_batch << endl;
    if (current_batch > 0) {
      edge_additions.del();
      edge_deletions.del();
      deletions_data.reset();
    }
  }

  void validateAndOpenFifo() {
    if (!stream_path) {
      std::cerr << "Error in command-line arguments. Please input a valid "
                   "stream_path or fileName"
                << std::endl;
      exit(1);
    }
    if (access(stream_path, F_OK) == -1) {
      int res = mkfifo(stream_path, 0777);
      if (res != 0) {
        std::cerr << "Could not create fifo" << std::endl;
        exit(1);
      }
    }
    cout << "Opening Stream: Waiting for writer to open..." << endl;
    stream_file.open(stream_path);
    cout << "Stream opened" << endl;
  }

  tuple<edgeArray, edgeArray, long, long>
  getNewEdgesFromFile(ifstream &inputFile, long numEdges, graph<vertex> GA,
                      bool symmetric, bool simpleFlag, bool fixedBatchFlag,
                      bool edgeValidityFlag, bool debugFlag,
                      bool &streamClosed) {
    if (numEdges == 0) {
#ifdef EDGEDATA
      return make_tuple(edgeArray(nullptr, nullptr, 0, 0),
                        edgeArray(nullptr, nullptr, 0, 0), 0, 0);
#else
      return make_tuple(edgeArray(nullptr, 0, 0), edgeArray(nullptr, 0, 0), 0,
                        0);
#endif
    }
    edge *EA = newA(edge, numEdges);
    edge *ED = newA(edge, numEdges);
    uintV source, destination;
#ifdef EDGEDATA
    EdgeData *edgeWeightEA = newA(EdgeData, numEdges);
    EdgeData *edgeWeightED = newA(EdgeData, numEdges);
    EdgeData *checkedEdgeWeightEA = newA(EdgeData, numEdges);
    EdgeData *checkedEdgeWeightED = newA(EdgeData, numEdges);
#endif
    char edgeType;
    string line;
    vector<string> tokens;
    long lineCount = 0;
    uintV maxVertex = 0;

#ifdef EDGEDATA
    intWeights *uncheckedEA = newA(intWeights, numEdges);
    intWeights *uncheckedED = newA(intWeights, numEdges);
#else
    intPair *uncheckedEA = newA(intPair, numEdges);
    intPair *uncheckedED = newA(intPair, numEdges);
#endif
    long uncheckedEACount = 0;
    long uncheckedEDCount = 0;
    long uncheckedEACountOrig = 0;
    long uncheckedEDCountOrig = 0;
    long edgesRead = numEdges;
    long edgesToRead = numEdges;

    long checkedEACount = 0;
    long checkedEDCount = 0;
    long cancelledEdges = 0;

    StreamEdge *edgesReceived = newA(StreamEdge, numEdges);
    cout << "Batch Size: " << numEdges << endl;
    do {
      edgesToRead = numEdges - checkedEDCount - checkedEACount;
      if (debugFlag) {
        cout << "Edges Added: " << checkedEACount << endl;
        cout << "Edges Deleted: " << checkedEDCount << endl;
        cout << "Edges to be Read: " << edgesToRead << endl;
      }
      uncheckedEACount = 0;
      uncheckedEDCount = 0;
      for (long i = 0; i < edgesToRead; i++) {
        long long numberOfBytesAvail = inputFile.rdbuf()->in_avail();
        if (!fixedBatchFlag && numberOfBytesAvail <= 0) {
          if (i == 0) {
            cout << "No Edges in Stream: Waiting for more edges or for stream "
                    "to close"
                 << endl;
          } else {
            edgesRead = i;
            cerr << "WARNING: Not enough edges to fulfill batch size. Only "
                 << i << " edges read." << endl;
            break;
          }
        }

        std::getline(inputFile, line);
        if ((line[0] == '%') || (line[0] == '#')) {
          continue;
        }

        if (!inputFile.good()) {
          edgesRead = i;
          streamClosed = true;
          cout << "WARNING: Stream Closed. Only " << i << " edges read" << endl;
          break;
        }
        tokens.clear();
        string buf;
        stringstream ss(line);
        while (ss >> buf) {
          tokens.push_back(buf);
        }
#ifdef EDGEDATA
        if (tokens.size() == 4) {
          edgeType = tokens[0].at(0);
          source = stoi(tokens[1]);
          destination = stoi(tokens[2]);

          if (edgeType == 'a') {
            new (edgeWeightEA + uncheckedEACount) EdgeData();
            edgeWeightEA[uncheckedEACount].createEdgeData(tokens[3].c_str());
            uncheckedEA[uncheckedEACount] =
                make_pair(source, make_pair(destination,
                                            &edgeWeightEA[uncheckedEACount]));
            uncheckedEACount++;
          } else if (edgeType == 'd') {
            new (edgeWeightED + uncheckedEDCount) EdgeData();
            edgeWeightED[uncheckedEDCount].createEdgeData(tokens[3].c_str());
            uncheckedED[uncheckedEDCount] =
                make_pair(source, make_pair(destination,
                                            &edgeWeightED[uncheckedEDCount]));
            uncheckedEDCount++;
          }
        }
#else
        if (tokens.size() == 3) {
          edgeType = tokens[0].at(0);
          source = stoi(tokens[1]);
          destination = stoi(tokens[2]);

          if (edgeType == 'a') {
            uncheckedEA[uncheckedEACount] = make_pair(source, destination);
            uncheckedEACount++;
          } else if (edgeType == 'd') {
            uncheckedED[uncheckedEDCount] = make_pair(source, destination);
            uncheckedEDCount++;
          }
        }
#endif
        else {
          std::cout << "Incorrect input format \n" << std::endl;
          inputFile.close();
          exit(1);
        }
      }
      uncheckedEACountOrig = uncheckedEACount;
      uncheckedEDCountOrig = uncheckedEDCount;

#ifdef EDGEDATA
      quickSort(uncheckedEA, uncheckedEACount, tripleBothCmp());
      quickSort(uncheckedED, uncheckedEDCount, tripleBothCmp());
#else
      quickSort(uncheckedEA, uncheckedEACount, pairBothCmp<uintE>());
      quickSort(uncheckedED, uncheckedEDCount, pairBothCmp<uintE>());
#endif

      bool *EAflag = newAWithZero(bool, uncheckedEACount);
      bool *EDflag = newAWithZero(bool, uncheckedEDCount);

      // remove duplicates
      if (simpleFlag) {
        uncheckedEACount = removeDuplicates(uncheckedEA, uncheckedEACount,
                                            symmetric, debugFlag);

        uncheckedEDCount = removeDuplicates(uncheckedED, uncheckedEDCount,
                                            symmetric, debugFlag);
      }

      if (!fixedBatchFlag) {
        long eaIndex = 0;
        long edIndex = 0;

        // remove values that cancel out
        while (eaIndex < uncheckedEACount && edIndex < uncheckedEDCount) {
#ifdef EDGEDATA
          if (uncheckedEA[eaIndex].first == uncheckedED[edIndex].first &&
              uncheckedEA[eaIndex].second.first ==
                  uncheckedED[edIndex].second.first)
#else
          if (uncheckedEA[eaIndex].first == uncheckedED[edIndex].first &&
              uncheckedEA[eaIndex].second == uncheckedED[edIndex].second)
#endif
          {
            EAflag[eaIndex] = true;
            EDflag[edIndex] = true;
            if (debugFlag) {
#ifdef EDGEDATA
              cerr << "CANCELLED: " << uncheckedED[edIndex].first << "\t"
                   << uncheckedED[edIndex].second.first << "\t"
                   << uncheckedED[edIndex].second.second << "\n";
#else
              cerr << "CANCELLED: " << uncheckedED[edIndex].first << "\t"
                   << uncheckedED[edIndex].second << "\n";
#endif
            }
            eaIndex++;
            edIndex++;
            cancelledEdges++;
          } else if (uncheckedEA[eaIndex].first == uncheckedED[edIndex].first) {
#ifdef EDGEDATA
            if (uncheckedEA[eaIndex].second.first <
                uncheckedED[edIndex].second.first)
#else
            if (uncheckedEA[eaIndex].second < uncheckedED[edIndex].second)
#endif
            {
              eaIndex++;
            } else {
              edIndex++;
            }
          } else if (uncheckedEA[eaIndex].first < uncheckedED[edIndex].first) {
            eaIndex++;
          } else {
            edIndex++;
          }
        }
      }

      // remove all EdgeDeletions >= G.n as they don't exist in the graph
      parallel_for(long i = 0; i < uncheckedEDCount; i++) {
#ifdef EDGEDATA
        if (uncheckedED[i].first >= GA.n || uncheckedED[i].second.first >= GA.n) {
#else
        if (uncheckedED[i].first >= GA.n || uncheckedED[i].second >= GA.n) {
#endif
          EDflag[i] = true;
        }
      }

      if (simpleFlag) {
        // check if edge additions edge is already in initial graph
        // we don't want the same edge from two vertices
        parallel_for(long i = 0; i < uncheckedEACount; i++) {
          if (EAflag[i] == false && uncheckedEA[i].first < GA.n) {
            vertex sourceV = GA.V[uncheckedEA[i].first];
            parallel_for(uintV j = 0; j < sourceV.getOutDegree(); j++) {
#ifdef EDGEDATA
              if (sourceV.getOutNeighbor(j) == uncheckedEA[i].second.first)
#else
              if (sourceV.getOutNeighbor(j) == uncheckedEA[i].second)
#endif
              {
                EAflag[i] = true;
                if (debugFlag) {
#ifdef EDGEDATA
                  cerr << "INVALID: " << uncheckedEA[i].first << "\t"
                       << uncheckedEA[i].second.first << "\t"
                       << uncheckedEA[i].second.second << "\n";
#else
                  cerr << "INVALID: " << uncheckedEA[i].first << "\t"
                       << uncheckedEA[i].second << "\n";
#endif
                }
              }
            }
          }
        }
      }

      if (edgeValidityFlag) {
        parallel_for(long i = 0; i < uncheckedEDCount; i++) {
          // check if edge deletions is valid
          if (EDflag[i] == false) {
            EDflag[i] = true;
            vertex sourceV = GA.V[uncheckedED[i].first];
            // check if this edge is present in the graph
            parallel_for(uintV j = 0; j < sourceV.getOutDegree(); j++) {
#ifdef EDGEDATA
              if (sourceV.getOutNeighbor(j) == uncheckedED[i].second.first)
#else
              if (sourceV.getOutNeighbor(j) == uncheckedED[i].second)
#endif
              {
                EDflag[i] = false;
              }
            }
            if (EDflag[i] == true) {
              if (debugFlag) {
#ifdef EDGEDATA
                cerr << "INVALID: " << uncheckedED[i].first << "\t"
                     << uncheckedED[i].second.first << "\t"
                     << uncheckedED[i].second.second << "\n";
#else
                cerr << "INVALID: " << uncheckedED[i].first << "\t"
                     << uncheckedED[i].second << "\n";
#endif
              }
            }
          }
        }
      }

      long maxCount = max(uncheckedEACount, uncheckedEDCount);
      for (long index = 0; index < maxCount; index++) {
        if (index < uncheckedEACount && !EAflag[index]) {
          EA[checkedEACount].source = uncheckedEA[index].first;
#ifdef EDGEDATA
          EA[checkedEACount].destination = uncheckedEA[index].second.first;
          new (checkedEdgeWeightEA + checkedEACount) EdgeData();
          EA[checkedEACount].edgeData = &checkedEdgeWeightEA[checkedEACount];
          EA[checkedEACount].edgeData->setEdgeDataFromPtr(
              uncheckedEA[index].second.second);
#else
          EA[checkedEACount].destination = uncheckedEA[index].second;
#endif
          if (EA[checkedEACount].source > maxVertex)
            maxVertex = EA[checkedEACount].source;
          if (EA[checkedEACount].destination > maxVertex)
            maxVertex = EA[checkedEACount].destination;
          checkedEACount++;
        }

        if (index < uncheckedEDCount && !EDflag[index]) {
          ED[checkedEDCount].source = uncheckedED[index].first;
#ifdef EDGEDATA
          ED[checkedEDCount].destination = uncheckedED[index].second.first;
          new (checkedEdgeWeightED + checkedEDCount) EdgeData();
          ED[checkedEDCount].edgeData = &checkedEdgeWeightED[checkedEDCount];
          ED[checkedEDCount].edgeData->setEdgeDataFromPtr(
              uncheckedED[index].second.second);
#else
          ED[checkedEDCount].destination = uncheckedED[index].second;
#endif
          checkedEDCount++;
        }
      }

      // ensure there are no duplicates within the edges to add/delete
      if (edgeValidityFlag && simpleFlag) {
        quickSort(EA, checkedEACount, edgeBothCmp());
        quickSort(ED, checkedEDCount, edgeBothCmp());
        checkedEACount = removeDuplicates(EA, checkedEACount, numEdges,
                                          symmetric, debugFlag);
        checkedEDCount = removeDuplicates(ED, checkedEDCount, numEdges,
                                          symmetric, debugFlag);
      }
      free(EAflag);
      free(EDflag);
#ifdef EDGEDATA
      parallel_for(long i = 0; i < uncheckedEDCountOrig; i++) {
        edgeWeightED[i].del();
      }
      parallel_for(long i = 0; i < uncheckedEACountOrig; i++) {
        edgeWeightEA[i].del();
      }
#endif
    } while (fixedBatchFlag && !streamClosed &&
             (checkedEACount + checkedEDCount) < numEdges);

    free(uncheckedEA);
    free(uncheckedED);
#ifdef EDGEDATA
    free(edgeWeightEA);
    free(edgeWeightED);
#endif
    free(edgesReceived);
#ifdef EDGEDATA
    return make_tuple(
        edgeArray(EA, checkedEdgeWeightEA, checkedEACount, maxVertex),
        edgeArray(ED, checkedEdgeWeightED, checkedEDCount, 0), edgesRead,
        cancelledEdges);
#else
    return make_tuple(edgeArray(EA, checkedEACount, maxVertex),
                      edgeArray(ED, checkedEDCount, 0), edgesRead,
                      cancelledEdges);
#endif
  }

  bool processNextBatch() {
    current_batch++;
    if (current_batch > number_of_batches) {
      cout << "Hit Max Batch Size" << endl;
      return false;
    }
    cleanup();

    timer timer1, timer2, fullTimer;
    fullTimer.start();
    parallel_for(uintV i = 0; i < n; i++) updated_vertices[i] = 0;

    long num_edges_read_from_file;
    long num_cancelled_edges;

    timer1.start();
    tie(edge_additions, edge_deletions_temp, num_edges_read_from_file,
        num_cancelled_edges) =
        getNewEdgesFromFile(stream_file, max_batch_size, my_graph,
                            my_graph.isSymmetric(), simple_flag,
                            fixed_batch_flag, enforce_edge_validity_flag,
                            debug_flag, stream_closed);
    cout << "Reading Time : " << timer1.stop() << endl;

    if (stream_closed && num_edges_read_from_file == 0) {
      cout << "No Edges in Batch" << endl;
      cout << "Ending Execution" << endl;
      return false;
    }

    double deletions_map_creation_time = 0;
    timer1.start();
    deletions_data.updateWithEdgesArray(edge_deletions_temp);
    deletions_map_creation_time = timer1.next();
    // cout << "deletions_map_creation_time : " << deletions_map_creation_time
    //      << "\n";
    edge_deletions =
        my_graph.deleteEdges(deletions_data, updated_vertices, debug_flag);
    edge_deletions_temp.del();
    cout << "Edge deletion time : "
         << deletions_map_creation_time + timer1.next() << "\n";

    timer1.start();
    if (edge_additions.maxVertex >= n) {
      long n_new = edge_additions.maxVertex + 1;
      updated_vertices = renewA(bool, updated_vertices, n_new);
      parallel_for(uintV i = my_graph.n; i < n_new; i++) {
        updated_vertices[i] = 1;
      }
      // update deletions_data
      deletions_data.updateNumVertices(n_new);
    }
    my_graph.addVertices(edge_additions.maxVertex);
    edge_additions = my_graph.addEdges(edge_additions, updated_vertices);
    cout << "Edge addition time : " << timer1.next() << "\n";
    if ((edge_additions.size > 0) || (edge_deletions.size > 0)) {
      return true;
    }
    cout << "No Edges in Stream" << endl;
    return false;
  }
};
#endif
