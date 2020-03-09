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

#ifdef EDGEDATA
// NOTE: The edge data type header file should then be included as the first header
// file at the top of the user program.
#include "SSSP_edgeData.h"
#endif

#include "../core/common/utils.h"
#include "../core/graphBolt/KickStarterEngine.h"
#include "../core/main.h"
#include <math.h>

#define MAX_DISTANCE 65535

// ======================================================================
// SSSPINFO
// ======================================================================
class SsspInfo {
public:
  uintV source_vertex;
  long weight_cap;

  SsspInfo() : source_vertex(0), weight_cap(3) {}

  SsspInfo(uintV _source_vertex, long _weight_cap)
      : source_vertex(_source_vertex), weight_cap(_weight_cap) {}

#ifdef EDGEDATA
#else
  uint16_t getWeight(uintV u, uintV v) {
    return (uint16_t)((u + v) % weight_cap + 1);
  }
#endif

  void copy(const SsspInfo &object) {
    source_vertex = object.source_vertex;
    weight_cap = object.weight_cap;
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {}

  void cleanup() {}
};

// ======================================================================
// VERTEXVALUE INITIALIZATION
// ======================================================================
template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  if (v != global_info.source_vertex) {
    v_vertex_value = MAX_DISTANCE;
  } else {
    v_vertex_value = 0;
  }
}

// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR FIRST ITERATION
// ======================================================================
template <class GlobalInfoType>
inline bool frontierVertex(const uintV &v, const GlobalInfoType &global_info) {
  if (v == global_info.source_vertex) {
    return true;
  } else {
    return false;
  }
}

// ======================================================================
// EDGE FUNCTION
// ======================================================================
template <class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline bool
edgeFunction(const uintV &u, const uintV &v, const EdgeDataType &edge_weight,
             const VertexValueType &u_value, VertexValueType &v_value,
             GlobalInfoType &global_info) {
  if (u_value == MAX_DISTANCE) {
    return false;
  } else {
#ifdef EDGEDATA
    v_value = u_value + edge_weight.weight;
#else
    v_value = u_value + global_info.getWeight(u, v);
#endif
    return true;
  }
}

// ======================================================================
// SHOULDPROPAGATE
// ======================================================================
// shouldPropagate condition for deciding if the value change in
// updated graph violates monotonicity
template <class VertexValueType, class GlobalInfoType>
inline bool shouldPropagate(const VertexValueType &old_value,
                            const VertexValueType &new_value,
                            GlobalInfoType &global_info) {
  return (new_value > old_value);
}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info) {}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  long n = G.n;
  int source_vertex = config.getOptionLongValue("-source", 0);
  int weight_cap = config.getOptionLongValue("-weight_cap", 5);
  SsspInfo global_info(source_vertex, weight_cap);

  cout << "Initializing engine ....\n";
  KickStarterEngine<vertex, uint16_t, SsspInfo> engine(G, global_info, config);
  engine.init();
  cout << "Finished initializing engine\n";
  engine.run();
}
