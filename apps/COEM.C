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
// NOTE: The edge data type header file should then be included as the first
// header file at the top of the user program.
#include "CF_edgeData.h"
#endif

#include "../core/common/matrix.h"
#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_simple.h"
#include "../core/main.h"
#include <math.h>

#define PARTITION_FILE_DEFAULT ""
#define SEED_FILE_DEFAULT ""

#define MODDER 0.5
#define CONST1 1.7777777777
#define CONST2 0.1

// ======================================================================
// CFINFO
// ======================================================================
template <class vertex> class CoemInfo {
public:
  // Should I just use vectors for this?
  uintV n;
  double epsilon;
  graph<vertex> *my_graph;

  double *in_weights;
  bool *partition_flags;
  bool *seed_flags;
  uintV flag_arrays_size;

  CoemInfo() : my_graph(nullptr), n(0), epsilon(0.01), flag_arrays_size(0) {}

  CoemInfo(graph<vertex> *_my_graph, uintV _n, double _epsilon)
      : my_graph(_my_graph), n(_n), epsilon(_epsilon) {
    if (n > 0) {
      flag_arrays_size = n;
      in_weights = newA(double, flag_arrays_size);
      partition_flags = newA(bool, flag_arrays_size);
      seed_flags = newA(bool, flag_arrays_size);
      parallel_for(uintV i = 0; i < n; i++) { in_weights[i] = 0.0; }
    }
  }

  void init() { computeInWeights(); }

#ifdef EDGEDATA
#else
  inline double getWeight(long i, long j) {
    return (fmod((i + j) * CONST1, MODDER) + CONST2);
  }
#endif

  // Helper functions for setting up partitions for non-bi-partite graphs
  void setPartitionsFromFile(string partitions_file_path) {
    if (partitions_file_path.compare(PARTITION_FILE_DEFAULT) == 0) {
      cout << "ERROR : INCORRECT Partitions FILE"
           << "\n";
      exit(1);
    }
    char temp[partitions_file_path.length() + 1];
    strcpy(temp, partitions_file_path.c_str());
    temp[partitions_file_path.length()] = '\0';
    _seq<char> S = readStringFromFile(temp);
    words W = stringToWords(S.A, S.n);
    long len = W.m;
    // Initialize to 0
    parallel_for(uintV i = 0; i < flag_arrays_size; i++) {
      partition_flags[i] = 0;
    }
    parallel_for(long i = 0; i < len; i++) {
      long vertexId = atoll(W.Strings[i]);
      if (vertexId > flag_arrays_size) {
        cout << "ERROR : " << vertexId << "\n";
      }
      partition_flags[vertexId] = 1;
    }
    W.del();
  }

  void setSeedsFromFile(string seeds_file_path) {
    if (seeds_file_path.compare(SEED_FILE_DEFAULT) == 0) {
      cout << "ERROR : INCORRECT Seeds FILE"
           << "\n";
      exit(1);
    }
    char temp[seeds_file_path.length() + 1];
    strcpy(temp, seeds_file_path.c_str());
    temp[seeds_file_path.length()] = '\0';
    _seq<char> S = readStringFromFile(temp);
    words W = stringToWords(S.A, S.n);
    long len = W.m;
    // Initialize to 0
    parallel_for(uintV i = 0; i < n; i++) { seed_flags[i] = 0; }
    parallel_for(long i = 0; i < len; i++) {
      long vertex_id = atol(W.Strings[i]);
      if (vertex_id > n) {
        cout << "ERROR : " << vertex_id << "\n";
      }
      if (belongsToNamesPartition(vertex_id)) {
        seed_flags[vertex_id] = 1;
      }
    }
    W.del();
  }

  inline bool isSeed(uintV i) const {
    return (i < flag_arrays_size) ? seed_flags[i] : 0;
  }

  inline bool belongsToNamesPartition(uintV i) const {
    return (i < flag_arrays_size) ? partition_flags[i] : i % 2;
  }

  inline bool belongsToCPartition(uintV i) const {
    return (i < flag_arrays_size) ? (!partition_flags[i]) : !(i % 2);
  }

  CoemInfo &operator=(const CoemInfo &object) {
    // copy other static objects
    n = object.n;
    epsilon = object.epsilon;
    flag_arrays_size = object.flag_arrays_size;
    in_weights = object.in_weights;
    partition_flags = object.partition_flags;
    seed_flags = object.seed_flags;
  }

  void copy(const CoemInfo &object) {
    // copy other static objects
    // We maintain a separate array for in_weights. The same array can be used
    // for the partition_flags and seed_flags.
    if (object.n > n) {
      if (n == 0) {
        n = object.n;
        in_weights = newA(double, n);
      } else {
        // realloc
        n = object.n;
        in_weights = renewA(double, in_weights, n);
      }
    }
    // copy the in_weights
    my_graph = object.my_graph;
    long min_n = std::min(object.n, n);
    parallel_for(uintV i = 0; i < n; i++) {
      in_weights[i] = object.in_weights[i];
    }

    epsilon = object.epsilon;
    flag_arrays_size = object.flag_arrays_size;
    partition_flags = object.partition_flags;
    seed_flags = object.seed_flags;
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    if (edge_additions.maxVertex >= n) {
      uintV n_old = n;
      n = edge_additions.maxVertex + 1;

      in_weights = renewA(double, in_weights, n);
      partition_flags = renewA(bool, partition_flags, n);
      seed_flags = renewA(bool, seed_flags, n);
      parallel_for(uintV i = n_old; i < n; i++) {
        partition_flags[i] = i % 2;
        seed_flags[i] = 0;
        in_weights[i] = 0;
      }
      flag_arrays_size = n;
    }
    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;
      double curr_weight = getWeight(source, destination);
      writeAdd(&in_weights[destination], curr_weight);
    }
    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;
      double curr_weight = getWeight(source, destination);
      writeAdd(&in_weights[destination], -curr_weight);
    }
  }

  void computeInWeights() {
    // cout << "computeInWeights\n";
    parallel_for(uintV v = 0; v < n; v++) {
      in_weights[v] = 0;
      intE inDegree = my_graph->V[v].getInDegree();
      for (long j = 0; j < inDegree; j++) {
        uintV u = my_graph->V[v].getInNeighbor(j);
        in_weights[v] += getWeight(u, v);
      }
    }
  }

  void cleanup() {
    if (flag_arrays_size > 0) {
      free(partition_flags);
      free(seed_flags);
    }
  }

  ~CoemInfo() {
    if (flag_arrays_size > 0) {
      free(in_weights);
    }
  }
};

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
double initial_aggregation_value = 0;
double initial_vertex_value = 0;
double initial_vertex_value_for_seeds = 1;
double aggregation_value_identity = 0;
double vertex_value_identity = 0;

template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const uintV &v,
                                       AggregationValueType &v_vertex_value,
                                       const GlobalInfoType &global_info) {
  v_vertex_value = initial_aggregation_value;
}

template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  if (global_info.isSeed(v)) {
    v_vertex_value = initial_vertex_value_for_seeds;
  } else {
    v_vertex_value = initial_vertex_value;
  }
}

template <class AggregationValueType>
inline AggregationValueType &aggregationValueIdentity() {
  return aggregation_value_identity;
}

template <class VertexValueType> inline VertexValueType &vertexValueIdentity() {
  return vertex_value_identity;
}

// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR A GIVEN ITERATION
// ======================================================================
template <class GlobalInfoType>
inline bool forceActivateVertexForIteration(const uintV &v, int iter,
                                            const GlobalInfoType &global_info) {
  bool ret = false;
  if ((iter == 1) && global_info.isSeed(v)) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}

template <class GlobalInfoType>
inline bool forceComputeVertexForIteration(const uintV &v, int iter,
                                           const GlobalInfoType &global_info) {
  return false;
}

inline bool shouldUseDelta(int iter) {
  // if ((iter == 1) || (iter == 2)) {
  if (iter == 1) {
    return false;
  }
  return true;
}

// ======================================================================
// ADD TO OR REMOVE FROM AGGREGATION VALUES
// ======================================================================
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregation(const AggregationValueType &incoming_value,
                             AggregationValueType &aggregate_value,
                             GlobalInfoType &global_info) {
  aggregate_value += incoming_value;
}
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, incoming_value);
}

template <class AggregationValueType, class GlobalInfoType>

inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info) {
  aggregate_value -= incoming_value;
}
template <class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info) {
  writeAdd(&aggregate_value, -incoming_value);
}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void computeFunction(const uintV &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_curr,
                            VertexValueType &vertex_value_next,
                            GlobalInfoType &global_info) {
  if (global_info.isSeed(v)) {
    vertex_value_next = vertex_value_curr;
    return;
  }

  if (global_info.my_graph->V[v].getInDegree() != 0) {
    vertex_value_next = aggregation_value / global_info.in_weights[v];
  } else {
    vertex_value_next = 0;
  }
  return;
}

template <class VertexValueType, class GlobalInfoType>
inline bool notDelZero(const VertexValueType &value_curr,
                      const VertexValueType &value_next,
                      GlobalInfoType &global_info) {
  if (fabs(value_next - value_curr) > global_info.epsilon) {
    return true;
  }
  return false;
}

// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class StaticInfoType>
inline void sourceChangeInContribution(
    const uintV &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev, const VertexValueType &v_value_curr,
    StaticInfoType &global_info) {
  v_change_in_contribution = v_value_curr - v_value_prev;
}

template <class AggregationValueType, class VertexValueType, class EdgeDataType,
          class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {
  if (global_info.belongsToNamesPartition(u) ==
      global_info.belongsToNamesPartition(v)) {
    return false;
  }
#ifdef EDGEDATA
  double edge_weight_val = edge_weight.weight;
#else
  double edge_weight_val = global_info.getWeight(u, v);
#endif

  u_change_in_contribution = edge_weight_val * u_change_in_contribution;

  return true;
}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// ======================================================================
template <class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     bool &forceComputeInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {
  if (global_info.in_weights[v] != global_info_old.in_weights[v])
    forceComputeInCurrentIteration = true;
}
template <class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const uintV &v,
                                          UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          bool &forceComputeInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old) {

  if (global_info.in_weights[v] != global_info_old.in_weights[v])
    forceComputeInCurrentIteration = true;
}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
void printHistory(const uintV &v, AggregationValueType **agg_values,
                  VertexValueType **actual_values, GlobalInfoType &info,
                  int history_iterations) {
  for (int iter = 0; iter < history_iterations; iter++) {
    cout << iter << ",";
    cout << agg_values[iter][v] << ",";
    cout << actual_values[iter][v] << ",";
    cout << "\n";
  }
}
template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const uintV &v,
                         GlobalInfoType &info) {
  output_file << info.isSeed(v) << " ";
}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  // cout << setprecision(TIME_PRECISION);
  uintV n = G.n;
  // cout << setprecision(VAL_PRECISION);
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  max_iters += 1;
  string partitions_file_path =
      config.getOptionValue("-partitionsFile", PARTITION_FILE_DEFAULT);
  string seeds_file_path =
      config.getOptionValue("-seedsFile", PARTITION_FILE_DEFAULT);

  double epsilon = 0.01d;

  CoemInfo<vertex> global_info(&G, n, epsilon);

  cout << "Reading partitions file ....\n";
  global_info.setPartitionsFromFile(partitions_file_path);

  cout << "Reading seeds file ....\n";
  global_info.setSeedsFromFile(seeds_file_path);
  // global_info.computeInWeights();

  cout << "Initializing engine ....\n";
  GraphBoltEngineSimple<vertex, double, double, CoemInfo<vertex>> engine(
      G, max_iters, global_info, false, config);
  engine.init();
  cout << "Finished init\n";
  engine.run();
}
