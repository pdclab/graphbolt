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
#include "LP_edgeData.h"
#endif

#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_simple.h"
#include "../core/main.h"
#include <math.h>

#define DEBUG_ROUND 100
#define PARTITION_FILE_DEFAULT ""
#define SEED_FILE_DEFAULT ""
#define NUMBER_OF_FEATURES 2
#define MOD_VAL_LP 0.5d

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE STRUCTURES
// ======================================================================
class LPVertexAggregationData {
public:
  double ngh_sum[NUMBER_OF_FEATURES];
  LPVertexAggregationData() {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      ngh_sum[i] = 0;
    }
  }
  LPVertexAggregationData &operator=(const LPVertexAggregationData &object) {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      ngh_sum[i] = object.ngh_sum[i];
    }
  }
  friend ostream &operator<<(ostream &os, const LPVertexAggregationData &dt);
};

ostream &operator<<(ostream &os,
                    const LPVertexAggregationData &aggregation_data) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    os << aggregation_data.ngh_sum[i];
    if (i + 1 != NUMBER_OF_FEATURES)
      os << " ";
  }
  return os;
}

class LPVertexData {
public:
  double features[NUMBER_OF_FEATURES];
  LPVertexData() {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      features[i] = 0;
    }
  }
  LPVertexData &operator=(const LPVertexData &object) {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      features[i] = object.features[i];
    }
  }
  friend ostream &operator<<(ostream &os, const LPVertexData &dt);
};

ostream &operator<<(ostream &os, const LPVertexData &vdata) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    os << vdata.features[i];
    if (i + 1 != NUMBER_OF_FEATURES)
      os << " ";
  }
  return os;
}

// ======================================================================
// LPINFO
// ======================================================================
template <class vertex> class LPInfo {
public:
  // Should I just use vectors for this?
  graph<vertex> *g;
  uintV n;
  double epsilon;
  double mod_val;
  double alpha;

  bool *seed_flags;

  LPInfo() : g(nullptr), n(0), epsilon(0), mod_val(0), alpha(0) {}

  LPInfo(graph<vertex> *_g, uintV _n, double _epsilon, double _mod_val,
         double _alpha)
      : g(_g), n(_n), epsilon(_epsilon), mod_val(_mod_val), alpha(_alpha) {
    if (n > 0) {
      seed_flags = newA(bool, n);
    }
  }

  LPInfo(LPInfo &object) {
    g = object.g;
    n = object.n;
    epsilon = object.epsilon;
    mod_val = object.mod_val;
    alpha = object.alpha;
    seed_flags = object.seed_flags;
  }

#ifdef EDGEDATA
#else
  inline double getWeight(uintV i, uintV j) {
    return fmod((i + j) * 1.7777777777, mod_val);
  }
#endif

  inline double getFeature(uintV v, int features) const {
    return fmod((v + features + 3) * 1.23456, mod_val);
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
      seed_flags[vertex_id] = 1;
    }
    W.del();
  }

  inline bool isSeed(uintEE i) const { return seed_flags[i]; }

  void copy(const LPInfo &object) {
    // copy other static objects
    n = object.n;
    g = object.g;
    epsilon = object.epsilon;
    alpha = object.alpha;
    mod_val = object.mod_val;
    seed_flags = object.seed_flags;
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    if (edge_additions.maxVertex >= n) {
      uintV n_old = n;
      n = edge_additions.maxVertex + 1;
      seed_flags = renewA(bool, seed_flags, n);
      parallel_for(uintV i = n_old; i < n; i++) { seed_flags[i] = 0; }
    }

    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;
      // Incrementally update static_data
    }
    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;
      // Incrementally update static_data
    }
  }

  void cleanup() {
    if (n > 0)
      deleteA(seed_flags);
  }

  ~LPInfo() {}
};

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
LPVertexAggregationData initial_aggregation_value;
LPVertexAggregationData aggregation_value_identity;
LPVertexData actual_value_identity;

template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const uintV &v,
                                       AggregationValueType &v_actual_value,
                                       const GlobalInfoType &global_info) {
  v_actual_value = initial_aggregation_value;
}

template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const uintV &v,
                                  VertexValueType &v_actual_value,
                                  const GlobalInfoType &global_info) {
  double sum = 0;
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    v_actual_value.features[i] = global_info.getFeature(v, i);
    sum += v_actual_value.features[i];
  }
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    v_actual_value.features[i] = v_actual_value.features[i] / sum;
  }
}

template <class AggregationValueType>
inline AggregationValueType &aggregationValueIdentity() {
  return aggregation_value_identity;
}

template <class VertexValueType> inline VertexValueType &vertexValueIdentity() {
  return actual_value_identity;
}

// ======================================================================
// ACTIVATE VERTEX/COMPUTE VERTEX FOR A GIVEN ITERATION
// ======================================================================
template <class GlobalInfoType>
inline bool forceActivateVertexForIteration(const uintV &v, int iter,
                                            const GlobalInfoType &global_info) {

  bool ret = false;
  if (iter == 1) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}

template <class GlobalInfoType>
inline bool forceComputeVertexForIteration(const uintV &v, int iter,
                                           const GlobalInfoType &global_info) {
  return true;
}

inline bool shouldUseDelta(int iter) {
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
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    aggregate_value.ngh_sum[i] += incoming_value.ngh_sum[i];
  }
}
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    writeAdd(&aggregate_value.ngh_sum[i], incoming_value.ngh_sum[i]);
  }
}

template <class AggregationValueType, class GlobalInfoType>

inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    aggregate_value.ngh_sum[i] -= incoming_value.ngh_sum[i];
  }
}
template <class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    writeAdd(&aggregate_value.ngh_sum[i], -incoming_value.ngh_sum[i]);
  }
}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void computeFunction(const uintV &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &actual_value_curr,
                            VertexValueType &actual_value_next,
                            GlobalInfoType &global_info) {

  if (global_info.isSeed(v)) {
    actual_value_next = actual_value_curr;
    return;
  }

  if (global_info.g->V[v].getInDegree() == 0) {
    actual_value_next = actual_value_curr;
    return;
  }

  double sum = 0;
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    sum += aggregation_value.ngh_sum[i];
  }

  double sum2 = 0;
  if (sum != 0) {

    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      actual_value_next.features[i] =
          actual_value_curr.features[i] * global_info.alpha +
          (aggregation_value.ngh_sum[i] / sum) * (1 - global_info.alpha);
      sum2 += actual_value_next.features[i];
    }
  } else {

    actual_value_next = actual_value_curr;
    return;
  }
  if (sum2 != 0) {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      actual_value_next.features[i] = actual_value_next.features[i] / sum2;
    }
  } else {
    for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
      actual_value_next.features[i] = actual_value_curr.features[i];
    }
  }
}

template <class VertexValueType, class GlobalInfoType>
inline bool isChanged(const VertexValueType &value_curr,
                      const VertexValueType &value_next,
                      GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    if (fabs(value_next.features[i] - value_curr.features[i]) >
        global_info.epsilon) {
      return true;
    }
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
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    v_change_in_contribution.ngh_sum[i] =
        (v_value_curr.features[i] - v_value_prev.features[i]);
  }
}

template <class AggregationValueType, class VertexValueType, class EdgeDataType,
          class GlobalInfoType>
inline bool edgeFunction(const uintV &u, const uintV &v,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {
#ifdef EDGEDATA
  double weight = edge_weight.weight;
#else
  double weight = global_info.getWeight(u, v);
#endif
  for (int i = 0; i < NUMBER_OF_FEATURES; i++) {
    u_change_in_contribution.ngh_sum[i] =
        weight * u_change_in_contribution.ngh_sum[i];
  }
  return true;
}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// ======================================================================
template <class GlobalInfoType>
inline void hasSourceChangedByUpdate(const uintV &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {}

template <class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const uintV &v,
                                          UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old) {}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
void printHistory(const uintV &v, AggregationValueType **agg_values,
                  VertexValueType **actual_values, GlobalInfoType &info,
                  int history_iterations) {
  for (int iter = 0; iter < history_iterations; iter++) {
    cout << iter << "," << agg_values[iter][v] << "," << actual_values[iter][v]
         << "\n";
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
  uintV n = G.n;
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  max_iters += 1;
  // NOTE : The seeds file is a mandatory argument. Each seed vertex is specified in a separate line in the file.  Refer "setSeedsFromFile()" fucntion for setting the seed vertices.
  string seeds_file_path =
      config.getOptionValue("-seedsFile", PARTITION_FILE_DEFAULT);
  double mod_val = config.getOptionDoubleValue("-modVal", MOD_VAL_LP);
  double epsilon = 0.010000000000000000000000000d;
  double alpha = 0.150000000000000000000000000d;

  LPInfo<vertex> global_info(&G, n, epsilon, mod_val, alpha);

  cout << "Reading seeds file ....\n";
  global_info.setSeedsFromFile(seeds_file_path);

  cout << "Initializing engine ....\n";
  GraphBoltEngineSimple<vertex, LPVertexAggregationData, LPVertexData, LPInfo<vertex>>
      engine(G, max_iters, global_info, false, config);
  engine.init();
  cout << "Finished initializing engine\n";
  engine.run();
}
