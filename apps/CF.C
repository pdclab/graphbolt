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

#include "../core/common/matrix.h"
#include "../core/common/utils.h"
#include "../core/graphBolt/GraphBoltEngine_complex.h"
#include "../core/main.h"
#include <math.h>

#define DEBUG_ROUND 100
#define PARTITION_FILE_DEFAULT ""
#define PARITION_INIT_3 1
#define NUMBER_OF_FACTORS 2
#define DEFAULT_SEED 5

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE STRUCTURES
// ======================================================================
class CFVertexAggregationData {
public:
  double inverse_component[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  double second_component[NUMBER_OF_FACTORS];
  CFVertexAggregationData() {
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      second_component[i] = 0;
      for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
        inverse_component[i * NUMBER_OF_FACTORS + j] = 0;
      }
    }
  }
  CFVertexAggregationData &operator=(const CFVertexAggregationData &object) {
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      second_component[i] = object.second_component[i];
      for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
        inverse_component[i * NUMBER_OF_FACTORS + j] =
            object.inverse_component[i * NUMBER_OF_FACTORS + j];
      }
    }
  }
  friend ostream &operator<<(ostream &os, const CFVertexAggregationData &dt);
};

ostream &operator<<(ostream &os,
                    const CFVertexAggregationData &aggregation_data) {
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
      os << aggregation_data.inverse_component[i * NUMBER_OF_FACTORS + j]
         << ",";
    }
  }
  return os;
}

class CFVertexData {
public:
  double latent_factors[NUMBER_OF_FACTORS];
  CFVertexData() {
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      latent_factors[i] = 0;
    }
  }
  CFVertexData &operator=(const CFVertexData &object) {
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      latent_factors[i] = object.latent_factors[i];
    }
  }
  friend ostream &operator<<(ostream &os, const CFVertexData &dt);
};

ostream &operator<<(ostream &os, const CFVertexData &vdata) {
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    os << vdata.latent_factors[i];
    if (i + 1 != NUMBER_OF_FACTORS) {
      os << " ";
    }
  }
  return os;
}

// ======================================================================
// CFINFO 
// ======================================================================
class CFGlobalInfo {
public:
  // Should I just use vectors for this?
  long n;
  double lamda_identity_matrix[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  double epsilon;
  double mod_val;
  bool use_random_init;
  long random_init_seed;

  bool *partition_flags;
  long partition_flags_array_size;

  CFGlobalInfo()
      : n(0), epsilon(0), random_init_seed(DEFAULT_SEED), mod_val(MOD_VAL),
        use_random_init(false) {}

  CFGlobalInfo(long _n, double _epsilon, double _mod_val)
      : n(_n), epsilon(_epsilon), random_init_seed(DEFAULT_SEED),
        mod_val(_mod_val), use_random_init(false) {
    if (n > 0) {
      partition_flags_array_size = n;
      partition_flags = newA(bool, partition_flags_array_size);
    }
  }

  CFGlobalInfo(CFGlobalInfo &object) {
    n = object.n;
    epsilon = object.epsilon;
    random_init_seed = object.random_init_seed;
    use_random_init = object.use_random_init;
  }

  void useRandomInit(long t_seed) {
    use_random_init = true;
    random_init_seed = t_seed;
  }

  inline double rating(uintEE i, uintEE j) { return fmod((i + j), mod_val); }

  void createPartition(string partition_file_path) {
    setPartitionsFromFile(partition_file_path);
  }

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
    parallel_for(long i = 0; i < partition_flags_array_size; i++) {
      partition_flags[i] = 0;
    }
    parallel_for(long i = 0; i < len; i++) {
      long vertexId = atoll(W.Strings[i]);
      if (vertexId > partition_flags_array_size) {
        cout << "ERROR : " << vertexId << "\n";
      }
      partition_flags[vertexId] = 1;
    }
    W.del();
  }

  void setLambda(double lambda) {
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
        if (i == j) {
          lamda_identity_matrix[i * NUMBER_OF_FACTORS + j] = lambda;
        } else {
          lamda_identity_matrix[i * NUMBER_OF_FACTORS + j] = 0;
        }
      }
    }
  }

  inline bool belongsToPartition1(uintEE i) const {
    return (i < partition_flags_array_size) ? partition_flags[i] : i % 2;
  }

  inline bool belongsToPartition2(uintEE i) const {
    return (i < partition_flags_array_size) ? (!partition_flags[i]) : !(i % 2);
  }

  CFGlobalInfo &operator=(const CFGlobalInfo &object) {
    // copy other static objects
    n = object.n;
    epsilon = object.epsilon;
    mod_val = object.mod_val;
    random_init_seed = object.random_init_seed;
    use_random_init = object.use_random_init;
    for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
      for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
        lamda_identity_matrix[i * NUMBER_OF_FACTORS + j] =
            object.lamda_identity_matrix[i * NUMBER_OF_FACTORS + j];
      }
    }
    partition_flags = object.partition_flags;
    partition_flags_array_size = object.partition_flags_array_size;
  }

  void copy(const CFGlobalInfo &object) {
    // copy other static objects
    n = object.n;
    epsilon = object.epsilon;
    mod_val = object.mod_val;
    setLambda(object.lamda_identity_matrix[0]);
    random_init_seed = object.random_init_seed;
    use_random_init = object.use_random_init;
    partition_flags = object.partition_flags;
    partition_flags_array_size = object.partition_flags_array_size;
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    if (edge_additions.maxVertex >= n) {
      long n_old = n;
      n = edge_additions.maxVertex + 1;
    }
  }

  void cleanup() {
    if (partition_flags_array_size > 0)
      free(partition_flags);
  }

  ~CFGlobalInfo() {}
};

// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
CFVertexAggregationData initial_aggregation_value;
CFVertexAggregationData aggregation_value_identity;
CFVertexData vertex_value_identity;

template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const long &v,
                                       AggregationValueType &v_vertex_value,
                                       const GlobalInfoType &global_info) {
  v_vertex_value = initial_aggregation_value;
}

template <class VertexValueType, class GlobalInfoType>
inline void initializeVertexValue(const long &v,
                                  VertexValueType &v_vertex_value,
                                  const GlobalInfoType &global_info) {
  for (long i = 0; i < NUMBER_OF_FACTORS; i++) {
    if (global_info.use_random_init) {
      double temp1 =
          (double)(global_info.random_init_seed +
                   hashInt((unsigned long)(v * NUMBER_OF_FACTORS + i)));
      v_vertex_value.latent_factors[i] = fmod(temp1, global_info.mod_val);
    } else {
      v_vertex_value.latent_factors[i] = global_info.mod_val;
    }
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
inline bool forceActivateVertexForIteration(const long &v, int iter,
                                            const GlobalInfoType &global_info) {
  bool ret = false;
  if ((iter == 1) && global_info.belongsToPartition2(v)) {
    ret = true;
  } else if ((iter == 2) && global_info.belongsToPartition1(v)) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}

template <class GlobalInfoType>
inline bool forceComputeVertexForIteration(const long &v, int iter,
                                           const GlobalInfoType &global_info) {
  bool ret = false;
  if ((iter == 1) && global_info.belongsToPartition1(v)) {
    ret = true;
  } else if ((iter == 2) && global_info.belongsToPartition2(v)) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}

inline bool shouldUseDelta(int iter) {
  if ((iter == 1) || (iter == 2)) {
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
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    aggregate_value.second_component[i] += incoming_value.second_component[i];
    for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
      aggregate_value.inverse_component[i * NUMBER_OF_FACTORS + j] +=
          incoming_value.inverse_component[i * NUMBER_OF_FACTORS + j];
    }
  }
}
template <class AggregationValueType, class GlobalInfoType>
inline void addToAggregationAtomic(const AggregationValueType &incoming_value,
                                   AggregationValueType &aggregate_value,
                                   GlobalInfoType &global_info) {}

template <class AggregationValueType, class GlobalInfoType>

inline void removeFromAggregation(const AggregationValueType &incoming_value,
                                  AggregationValueType &aggregate_value,
                                  GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    aggregate_value.second_component[i] -= incoming_value.second_component[i];
    for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
      aggregate_value.inverse_component[i * NUMBER_OF_FACTORS + j] -=
          incoming_value.inverse_component[i * NUMBER_OF_FACTORS + j];
    }
  }
}
template <class AggregationValueType, class GlobalInfoType>
inline void
removeFromAggregationAtomic(const AggregationValueType &incoming_value,
                            AggregationValueType &aggregate_value,
                            GlobalInfoType &global_info) {}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline void computeFunction(const long &v,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_curr,
                            VertexValueType &vertex_value_next,
                            GlobalInfoType &global_info) {

  double inverse_component_with_lamda[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  double inverse_component_final[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
      inverse_component_with_lamda[i * NUMBER_OF_FACTORS + j] =
          aggregation_value.inverse_component[i * NUMBER_OF_FACTORS + j] +
          global_info.lamda_identity_matrix[i * NUMBER_OF_FACTORS + j];
    }
  }
  inverseOfMatrix(inverse_component_with_lamda, NUMBER_OF_FACTORS,
                  NUMBER_OF_FACTORS, inverse_component_final, NUMBER_OF_FACTORS,
                  NUMBER_OF_FACTORS);
  multiplyMatrices(inverse_component_final, NUMBER_OF_FACTORS,
                   NUMBER_OF_FACTORS, aggregation_value.second_component,
                   NUMBER_OF_FACTORS, 1, vertex_value_next.latent_factors,
                   NUMBER_OF_FACTORS, 1);
}

template <class VertexValueType, class GlobalInfoType>
inline bool isChanged(const VertexValueType &value_curr,
                         const VertexValueType &value_next,
                         GlobalInfoType &global_info) {
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    if (fabs(value_next.latent_factors[i] - value_curr.latent_factors[i]) >
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
    const long &v, AggregationValueType &v_change_in_contribution,
    const VertexValueType &v_value_prev, const VertexValueType &v_value_curr,
    StaticInfoType &global_info) {}

template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline bool edgeFunction(const long &u, const long &v,
                         const VertexValueType &u_value,
                         AggregationValueType &u_change_in_contribution,
                         GlobalInfoType &global_info) {
  if (global_info.belongsToPartition1(u) ==
      global_info.belongsToPartition1(v)) {
    return false;
  }
  double edge_rating = global_info.rating(u, v);
  double u_transpose_product[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  getTransposeProduct<double>(u_value.latent_factors, NUMBER_OF_FACTORS, 1,
                              u_change_in_contribution.inverse_component,
                              NUMBER_OF_FACTORS, NUMBER_OF_FACTORS);
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    u_change_in_contribution.second_component[i] =
        edge_rating * u_value.latent_factors[i];
  }
  return true;
}

template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
inline bool edgeFunctionDelta(const long &u, const long &v,
                              const VertexValueType &u_value_prev,
                              const VertexValueType &u_value_curr,
                              AggregationValueType &u_change_in_contribution,
                              GlobalInfoType &global_info) {
  if (global_info.belongsToPartition1(u) ==
      global_info.belongsToPartition1(v)) {
    return false;
  }
  double edge_rating = global_info.rating(u, v);
  double u_transpose_product_prev[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  double u_transpose_product_curr[NUMBER_OF_FACTORS * NUMBER_OF_FACTORS];
  getTransposeProduct<double>(u_value_curr.latent_factors, NUMBER_OF_FACTORS, 1,
                              u_transpose_product_curr, NUMBER_OF_FACTORS,
                              NUMBER_OF_FACTORS);
  getTransposeProduct<double>(u_value_prev.latent_factors, NUMBER_OF_FACTORS, 1,
                              u_transpose_product_prev, NUMBER_OF_FACTORS,
                              NUMBER_OF_FACTORS);
  for (int i = 0; i < NUMBER_OF_FACTORS; i++) {
    u_change_in_contribution.second_component[i] =
        edge_rating *
        (u_value_curr.latent_factors[i] - u_value_prev.latent_factors[i]);
    for (int j = 0; j < NUMBER_OF_FACTORS; j++) {
      u_change_in_contribution.inverse_component[i * NUMBER_OF_FACTORS + j] =
          (u_transpose_product_curr[i * NUMBER_OF_FACTORS + j] -
           u_transpose_product_prev[i * NUMBER_OF_FACTORS + j]);
    }
  }
  return true;
}

// ======================================================================
// INCREMENTAL COMPUTING / DETERMINING FRONTIER
// ======================================================================
template <class GlobalInfoType>
inline void hasSourceChangedByUpdate(const long &v, UpdateType update_type,
                                     bool &activateInCurrentIteration,
                                     GlobalInfoType &global_info,
                                     GlobalInfoType &global_info_old) {}
template <class GlobalInfoType>
inline void hasDestinationChangedByUpdate(const long &v, UpdateType update_type,
                                          bool &activateInCurrentIteration,
                                          GlobalInfoType &global_info,
                                          GlobalInfoType &global_info_old) {}

// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
void printHistory(const long &v, AggregationValueType **agg_values,
                  VertexValueType **actual_values, GlobalInfoType &info,
                  int history_iterations) {
  for (int iter = 0; iter < history_iterations; iter++) {
    cout << iter << ",";
    for (long i = 0; i < NUMBER_OF_FACTORS; i++) {
      for (long j = 0; j < NUMBER_OF_FACTORS; j++) {
        cout << agg_values[iter][v].inverse_component[i * NUMBER_OF_FACTORS + j]
             << ",";
      }
    }
    for (long i = 0; i < NUMBER_OF_FACTORS; i++) {
      cout << actual_values[iter][v].latent_factors[i] << ",";
    }
    cout << "\n";
  }
}
template <class GlobalInfoType>
void printAdditionalData(ofstream &output_file, const long &v,
                         GlobalInfoType &info) {
  output_file << info.belongsToPartition1(v) << " ";
}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  // cout << setprecision(TIME_PRECISION);
  long n = G.n;
  // cout << setprecision(VAL_PRECISION);
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  max_iters += 1;
  string partitions_file_path =
      config.getOptionValue("-partitionsFile", PARTITION_FILE_DEFAULT);
  double mod_val = config.getOptionDoubleValue("-modVal", MOD_VAL);
  bool rand_init = config.getOption("-randInit");
  double lambda = config.getOptionDoubleValue("-lambda", 0.0001); // lambda
  double epsilon = 0.010000000000000000000000000d;

  CFGlobalInfo global_info(n, epsilon, mod_val);

  global_info.useRandomInit(5);
  global_info.setLambda(lambda);
  cout << "Reading partitions file ....\n";
  global_info.createPartition(partitions_file_path);
  cout << "Initializing engine ....\n";
  GraphBoltEngineComplex<vertex, CFVertexAggregationData, CFVertexData,
                         CFGlobalInfo>
      engine(G, max_iters, global_info, true, config);
  engine.init();
  cout << "Finished init\n";
  engine.run();
}
