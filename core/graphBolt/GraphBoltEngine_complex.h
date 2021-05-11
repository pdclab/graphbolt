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

#ifndef GRAPHBOLT_ENGINE_COMPLEX_H
#define GRAPHBOLT_ENGINE_COMPLEX_H

#include "GraphBoltEngine.h"

// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
// Additional function for Complex engine.
// For complex aggregations, the contribution for the previous iteration and
// contribution for current iteration for that edge have to computed
// individually and only then its change_in_contribution can be updated. Note
// that in case of algorithms with multiple aggregation values per vertex, some
// of them may be simple aggregations. So, even for complex aggregations, the
// sourceChangeInContribution() function will be called before edgeFunction() /
// edgeFunctionDelta() to take advantage of such cases. For a vertex u, given
// its value in current iteraion and previous iteraion, and its computed
// change_in_contribution, update the change_in_contribution for a given edge
// (u, v). Return false if the value from u should not be included in the
// aggregation value of v. Return true otherwise. NOTE: The changes to
// change_in_contribution should be made in place. No LOCKS/CAS required for
// updating u_change_in_contribution
template <class AggregationValueType, class VertexValueType, class EdgeDataType,
          class GlobalInfoType>
inline bool edgeFunctionDelta(const uintV &u, const uintV &v,
                              const EdgeDataType &edge_weight,
                              const VertexValueType &u_value_prev,
                              const VertexValueType &u_value_curr,
                              AggregationValueType &u_change_in_contribution,
                              GlobalInfoType &global_info);

// ======================================================================
// GRAPHBOLTENGINECOMPLEX
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class GraphBoltEngineComplex
    : public GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                             GlobalInfoType> {
public:
  // use_source_contribution is false by default for complex aggregations. In
  // certain types of complex aggregations, some part of the computation is
  // common for a given source vertex in an edgeMap. In these cases, set this
  // flag to true and this optimization can be taken advantage of.
  AggregationValueType *source_change_in_contribution_old;
  GraphBoltEngineComplex(graph<vertex> &_my_graph, int _max_iter,
                         GlobalInfoType &_static_data, bool _use_lock,
                         commandLine _config)
      : GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>(_my_graph, _max_iter, _static_data,
                                        _use_lock, _config) {
    use_source_contribution = false;
  }

  // ======================================================================
  // TEMPORARY STRUCTURES USED BY THE COMPLEX ENGINE
  // ======================================================================
  void createTemporaryStructures() {
    if (use_source_contribution) {
      source_change_in_contribution_old = newA(AggregationValueType, n);
    }
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::createTemporaryStructures();
  }
  void resizeTemporaryStructures() {
    if (use_source_contribution) {
      source_change_in_contribution_old =
          renewA(AggregationValueType, source_change_in_contribution_old, n);
    }
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::resizeTemporaryStructures();
    initTemporaryStructures(n_old, n);
  }
  void freeTemporaryStructures() {
    if (use_source_contribution) {
      deleteA(source_change_in_contribution_old);
    }
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::freeTemporaryStructures();
  }
  void initTemporaryStructures() { initTemporaryStructures(0, n); }
  void initTemporaryStructures(long start_index, long end_index) {
    if (use_source_contribution) {
      parallel_for(long v = start_index; v < end_index; v++) {
        source_change_in_contribution_old[v] =
            aggregationValueIdentity<AggregationValueType>();
      }
    }
    GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                    GlobalInfoType>::initTemporaryStructures(start_index,
                                                             end_index);
  }

  // ======================================================================
  // TRADITIONAL INCREMENTAL COMPUTATION
  // ======================================================================
  int traditionalIncrementalComputation(int start_iteration) {
    timer iteration_timer, phase_timer;
    double misc_time, copy_time, phase_time, iteration_time;

    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool use_delta = true;
    if (frontier_curr_vs.numNonzeros() == 0) {
      converged_iteration = start_iteration;
    } else {
      for (int iter = start_iteration; iter < max_iterations; iter++) {
        // initialize timers
        {
          phase_timer.start();
          misc_time = 0;
          copy_time = 0;
        }

        // ========== COPY - Prepare curr iteration ==========
        if (iter > 0) {
          // Copy the aggregate and actual value from iter-1 to iter
          parallel_for(uintV v = 0; v < n; v++) {
            vertex_values[iter][v] = vertex_values[iter - 1][v];
            aggregation_values[iter][v] = aggregation_values[iter - 1][v];
            delta[v] = aggregationValueIdentity<AggregationValueType>();
          }
        }
        use_delta = shouldUseDelta(iter);
        phase_time = phase_timer.next();
        adaptive_executor.updateCopyTime(iter, phase_time);

        // ========== MISC - count active edges for AE ==========
        adaptive_executor.updateEdgesProcessed(iter, my_graph,
                                               frontier_curr_vs);
        phase_time = phase_timer.next();
        misc_time = phase_time;

        adaptive_executor.updateMiscTime(iter, phase_time);

        // ========== EDGE COMPUTATION ==========
        if ((use_source_contribution) && (iter == 1)) {
          // Compute source contribution for first iteration
          parallel_for(uintV u = 0; u < n; u++) {
            if (frontier_curr[u]) {
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  u, source_change_in_contribution[u],
                  vertexValueIdentity<VertexValueType>(),
                  vertex_values[iter - 1][u], global_info);
            }
          }
        }

        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // check for propagate and retract for the vertices.
            intE outDegree = my_graph.V[u].getOutDegree();
            granular_for(j, 0, outDegree, (outDegree > 1024), {
              uintV v = my_graph.V[u].getOutNeighbor(j);
              AggregationValueType contrib_change =
                  use_source_contribution
                      ? source_change_in_contribution[u]
                      : aggregationValueIdentity<AggregationValueType>();
              bool ret = false;
#ifdef EDGEDATA
              EdgeData *edge_data = my_graph.V[u].getOutEdgeData(j);
#else
              EdgeData *edge_data = &emptyEdgeData;
#endif
              if (use_delta) {
                // For first iteration, usually noDelta
                ret = edgeFunctionDelta(
                    u, v, *edge_data, vertex_values[iter - 2][u],
                    vertex_values[iter - 1][u], contrib_change, global_info);
              } else {
                ret = edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                                   contrib_change, global_info);
              }

              // lock if needed
              if (ret) {
                if (use_lock) {
                  vertex_locks[v].writeLock();
                  addToAggregation(contrib_change, delta[v], global_info);
                  vertex_locks[v].unlock();
                } else {
                  addToAggregationAtomic(contrib_change, delta[v], global_info);
                }
                if (!frontier_next[v])
                  frontier_next[v] = 1;
              }
            });
          }
        }

        phase_time = phase_timer.next();
        adaptive_executor.updateEdgeMapTime(iter, phase_time);

        // ========== VERTEX COMPUTATION ==========
        parallel_for(uintV v = 0; v < n; v++) {
          // Reset frontier for next iteration
          frontier_curr[v] = 0;
          if (frontier_next[v] ||
              forceComputeVertexForIteration(v, iter, global_info)) {
            frontier_next[v] = 0;

            // Update aggregation value and reset change received[v] (i.e.
            // delta[v])
            addToAggregation(delta[v], aggregation_values[iter][v],
                             global_info);
            delta[v] = aggregationValueIdentity<AggregationValueType>();

            // Calculate new_value based on the updated aggregation value
            VertexValueType new_value;
            computeFunction(v, aggregation_values[iter][v],
                            vertex_values[iter - 1][v], new_value, global_info);

            // Check if change is significant
            if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
              // change is significant. Update vertex_values
              vertex_values[iter][v] = new_value;
              // Set active for next iteration.
              frontier_curr[v] = 1;
            } else {
              // change is not significant. Copy vertex_values[iter-1]
              vertex_values[iter][v] = vertex_values[iter - 1][v];
            }
          }
          frontier_curr[v] =
              frontier_curr[v] ||
              forceActivateVertexForIteration(v, iter + 1, global_info);
          if (frontier_curr[v]) {
            if (use_source_contribution) {
              // update source_contrib for next iteration
              sourceChangeInContribution<AggregationValueType, VertexValueType,
                                         GlobalInfoType>(
                  v, source_change_in_contribution[v],
                  vertex_values[iter - 1][v], vertex_values[iter][v],
                  global_info);
            }
          }
        }
        phase_time = phase_timer.stop();
        adaptive_executor.updateVertexMapTime(iter, phase_time);

        vertexSubset temp_vs(n, frontier_curr);
        frontier_curr_vs = temp_vs;
        misc_time += phase_timer.next();
        iteration_time = iteration_timer.stop();

        if (ae_enabled && iter == 1) {
          adaptive_executor.setApproximateTimeForCurrIter(iteration_time);
        }

        // Convergence check
        converged_iteration = iter;
        if (frontier_curr_vs.isEmpty()) {
          break;
        }
      }
    }
    adaptive_executor.updateEquation(converged_iteration);
    return converged_iteration;
  }

  // ======================================================================
  // DELTACOMPUTE
  // ======================================================================
  void deltaCompute(edgeArray &edge_additions, edgeArray &edge_deletions) {
    timer iteration_timer, phase_timer, full_timer, pre_compute_timer;
    double misc_time, copy_time, phase_time, iteration_time, pre_compute_time;
    full_timer.start();

    // Handle newly added vertices
    n_old = n;
    if (edge_additions.maxVertex >= n) {
      processVertexAddition(edge_additions.maxVertex);
    }

    // Reset values before incremental computation
    parallel_for(uintV v = 0; v < n; v++) {
      frontier_curr[v] = 0;
      frontier_next[v] = 0;
      changed[v] = 0;
      retract[v] = 0;
      propagate[v] = 0;

      vertex_value_old_prev[v] = vertexValueIdentity<VertexValueType>();
      vertex_value_old_curr[v] = vertexValueIdentity<VertexValueType>();
      initializeVertexValue<VertexValueType>(v, vertex_value_old_next[v],
                                             global_info);

      delta[v] = aggregationValueIdentity<AggregationValueType>();
      if (use_source_contribution) {
        source_change_in_contribution[v] =
            aggregationValueIdentity<AggregationValueType>();
        source_change_in_contribution_old[v] =
            aggregationValueIdentity<AggregationValueType>();
      }
    }

    // ==================== UPDATE GLOBALINFO ===============================

    global_info_old.copy(global_info);

    global_info.processUpdates(edge_additions, edge_deletions);

    // ========== EDGE COMPUTATION - DIRECT CHANGES - for first iter ==========
    pre_compute_timer.start();
    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;

      hasSourceChangedByUpdate(source, edge_addition_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_addition_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {
        // Update frontier and changed values

        if (frontier_curr[source]) {
          retract[source] = frontier_curr[source];
          propagate[source] = frontier_curr[source];
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          retract[destination] = frontier_curr[destination];
          propagate[destination] = frontier_curr[destination];
          changed[destination] = true;
        }

        AggregationValueType contrib_change =
            aggregationValueIdentity<AggregationValueType>();
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
        }

// Do repropagate for edge source->destination.
#ifdef EDGEDATA
        EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret =
            edgeFunction(source, destination, *edge_data,
                         vertex_values[0][source], contrib_change, global_info);
        if (ret) {
          if (use_lock) {
            vertex_locks[destination].writeLock();
            addToAggregation(contrib_change, delta[destination],
                             global_info_old);
            vertex_locks[destination].unlock();
          } else {
            addToAggregationAtomic(contrib_change, delta[destination],
                                   global_info_old);
          }
          if (!changed[destination])
            changed[destination] = true;
        }
      }
    }

    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;

      hasSourceChangedByUpdate(source, edge_deletion_enum,
                               frontier_curr[source], changed[source],
                               global_info, global_info_old);
      hasSourceChangedByUpdate(destination, edge_deletion_enum,
                               frontier_curr[destination], changed[destination],
                               global_info, global_info_old);
      if (forceActivateVertexForIteration(source, 1, global_info_old)) {
        // Update frontier and changed values
        if (frontier_curr[source]) {
          retract[source] = frontier_curr[source];
          propagate[source] = frontier_curr[source];
          changed[source] = true;
        }
        if (frontier_curr[destination]) {
          retract[destination] = frontier_curr[destination];
          propagate[destination] = frontier_curr[destination];
          changed[destination] = true;
        }

        AggregationValueType contrib_change =
            aggregationValueIdentity<AggregationValueType>();
        if (use_source_contribution) {
          sourceChangeInContribution<AggregationValueType, VertexValueType,
                                     GlobalInfoType>(
              source, contrib_change, vertexValueIdentity<VertexValueType>(),
              vertex_values[0][source], global_info_old);
        }

// Do retract for edge source->destination
#ifdef EDGEDATA
        EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
        EdgeData *edge_data = &emptyEdgeData;
#endif
        bool ret =
            edgeFunction(source, destination, *edge_data,
                         vertex_values[0][source], contrib_change, global_info);
        if (ret) {
          if (use_lock) {
            vertex_locks[destination].writeLock();
            removeFromAggregation(contrib_change, delta[destination],
                                  global_info_old);
            vertex_locks[destination].unlock();
          } else {
            removeFromAggregationAtomic(contrib_change, delta[destination],
                                        global_info_old);
          }
          if (!changed[destination])
            changed[destination] = true;
        }
      }
    }
    pre_compute_time = pre_compute_timer.stop();

    // =============== INCREMENTAL COMPUTE - REFINEMENT START ================
    vertexSubset frontier_curr_vs(n, frontier_curr);
    bool should_switch_now = false;
    bool use_delta = true;

    if (ae_enabled && shouldSwitch(0, 0)) {
      should_switch_now = true;
    }
    for (int iter = 1; iter < max_iterations; iter++) {

      // Perform switch
      if (ae_enabled && should_switch_now) {
        converged_iteration = performSwitch(iter);
        break;
      }

      // initialize timers
      {
        iteration_timer.start();
        phase_timer.start();
        misc_time = 0;
        copy_time = 0;
      }

      // iteration_timer.start();
      use_delta = shouldUseDelta(iter);

      // ==================== COPY - PREPARE CURRENT ITERATION ====================
      {
        VertexValueType *temp1 = vertex_value_old_prev;
        vertex_value_old_prev = vertex_value_old_curr;
        vertex_value_old_curr = vertex_value_old_next;
        vertex_value_old_next = temp1;

        if (iter <= converged_iteration) {
          parallel_for(uintV v = 0; v < n; v++) {
            vertex_value_old_next[v] = vertex_values[iter][v];
          }
        } else {
          converged_iteration = performSwitch(iter);
          break;
        }
      }
      copy_time += phase_timer.next();

      // ========== EDGEMAP - TRANSITIVE CHANGES ==========
      if ((use_source_contribution) && (iter == 1)) {
        // Compute source contribution for first iteration
        parallel_for(uintV u = 0; u < n; u++) {
          if (frontier_curr[u]) {
            // compute source change in contribution
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, source_change_in_contribution[u],
                vertexValueIdentity<VertexValueType>(),
                vertex_values[iter - 1][u], global_info);

            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                u, source_change_in_contribution_old[u],
                vertexValueIdentity<VertexValueType>(),
                vertex_value_old_curr[u], global_info_old);
          }
        }
      }

      parallel_for(uintV u = 0; u < n; u++) {
        if (frontier_curr[u]) {
          // check for propagate and retract for the vertices.
          intE outDegree = my_graph.V[u].getOutDegree();
          granular_for(i, 0, outDegree, (outDegree > 1024), {
            uintV v = my_graph.V[u].getOutNeighbor(i);
            bool ret_old = false;
            bool ret = false;

            AggregationValueType to_retract =
                use_source_contribution
                    ? source_change_in_contribution_old[u]
                    : aggregationValueIdentity<AggregationValueType>();
            AggregationValueType to_propagate =
                use_source_contribution
                    ? source_change_in_contribution[u]
                    : aggregationValueIdentity<AggregationValueType>();

            // retract
            if (retract[u]) {
#ifdef EDGEDATA
              EdgeData *edge_data = my_graph.V[u].getOutEdgeData(i);
#else
              EdgeData *edge_data = &emptyEdgeData;
#endif
              if (use_delta) {
                ret_old = edgeFunctionDelta(
                    u, v, *edge_data, vertex_value_old_prev[u],
                    vertex_value_old_curr[u], to_retract, global_info_old);
              } else {
                // For first iteration, noDelta
                ret_old =
                    edgeFunction(u, v, *edge_data, vertex_value_old_curr[u],
                                 to_retract, global_info_old);
              }
            }

            // propagate
            if (propagate[u]) {
#ifdef EDGEDATA
              EdgeData *edge_data = my_graph.V[u].getOutEdgeData(i);
#else
              EdgeData *edge_data = &emptyEdgeData;
#endif
              if (use_delta) {
                ret = edgeFunctionDelta(
                    u, v, *edge_data, vertex_values[iter - 2][u],
                    vertex_values[iter - 1][u], to_propagate, global_info);
              } else {
                // For first iteration, noDelta
                ret = edgeFunction(u, v, *edge_data, vertex_values[iter - 1][u],
                                   to_propagate, global_info);
              }
            }

            if (ret || ret_old) {
              if (use_lock) {
                vertex_locks[v].writeLock();
                if (ret_old) {
                  removeFromAggregation(to_retract, delta[v], global_info_old);
                }
                if (ret) {
                  addToAggregation(to_propagate, delta[v], global_info);
                }
                vertex_locks[v].unlock();

              } else {
                removeFromAggregationAtomic(to_retract, delta[v],
                                            global_info_old);
                addToAggregationAtomic(to_propagate, delta[v], global_info);
              }

              if (!changed[v])
                changed[v] = 1;
            }
          });
        }
      }
      phase_time = phase_timer.next();

      // ========== VERTEX COMPUTATION  ==========
      bool use_delta_next_iteration = shouldUseDelta(iter + 1);
      parallel_for(uintV v = 0; v < n; v++) {
        if ((v >= n_old) && (changed[v] == false)) {
          changed[v] = forceComputeVertexForIteration(v, iter, global_info);
        }
        // changed vertices need to be processed
        if (changed[v]) {
          frontier_curr[v] = 0;
          retract[v] = 0;
          propagate[v] = 0;

          // delta has the current cumulative change for the vertex.
          // Update the aggregation value in history
          addToAggregation(delta[v], aggregation_values[iter][v], global_info);

          VertexValueType new_value;
          computeFunction(v, aggregation_values[iter][v],
                          vertex_values[iter - 1][v], new_value, global_info);

          if (forceActivateVertexForIteration(v, iter + 1, global_info)) {
            frontier_curr[v] = 1;
            propagate[v] = 1;
            retract[v] = 1;
          }

          // Update v_change based on updated graph
          if (notDelZero(new_value, vertex_values[iter - 1][v], global_info)) {
            // change is significant. Update vertex_values
            vertex_values[iter][v] = new_value;
            frontier_curr[v] = 1;
            propagate[v] = 1;
            if (use_source_contribution) {
              if (use_delta_next_iteration) {
                sourceChangeInContribution<AggregationValueType,
                                           VertexValueType, GlobalInfoType>(
                    v, source_change_in_contribution[v],
                    vertex_values[iter - 1][v], vertex_values[iter][v],
                    global_info);
              } else {
                sourceChangeInContribution<AggregationValueType,
                                           VertexValueType, GlobalInfoType>(
                    v, source_change_in_contribution[v],
                    vertexValueIdentity<VertexValueType>(),
                    vertex_values[iter][v], global_info);
              }
            }
          } else {
            // change is not significant. Copy vertex_values[iter-1]
            vertex_values[iter][v] = vertex_values[iter - 1][v];
          }

          // Update v_change based on initial graph
          if (notDelZero(vertex_value_old_next[v], vertex_value_old_curr[v],
                        global_info_old)) {
            // change is significant. Update v_change
            frontier_curr[v] = 1;
            retract[v] = 1;
            if (use_source_contribution) {
              if (use_delta_next_iteration) {
                sourceChangeInContribution<AggregationValueType,
                                           VertexValueType, GlobalInfoType>(
                    v, source_change_in_contribution_old[v],
                    vertex_value_old_curr[v], vertex_value_old_next[v],
                    global_info_old);
              } else {

                sourceChangeInContribution<AggregationValueType,
                                           VertexValueType, GlobalInfoType>(
                    v, source_change_in_contribution_old[v],
                    vertexValueIdentity<VertexValueType>(),
                    vertex_value_old_next[v], global_info_old);
              }
            }
          }
        }
      }
      phase_time = phase_timer.next();

      // ========== EDGE COMPUTATION - DIRECT CHANGES - for next iter ==========
      bool has_direct_changes = false;
      parallel_for(long i = 0; i < edge_additions.size; i++) {
        uintV source = edge_additions.E[i].source;
        uintV destination = edge_additions.E[i].destination;

        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
          AggregationValueType contrib_change_old =
              aggregationValueIdentity<AggregationValueType>();

          // Do repropagate for edge source->destination.
          bool ret = false;
          if (use_delta_next_iteration) {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change_old, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
#ifdef EDGEDATA
            EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunctionDelta(source, destination, *edge_data,
                                    vertex_value_old_curr[source],
                                    vertex_value_old_next[source],
                                    contrib_change_old, global_info_old);
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change_old,
                vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
#ifdef EDGEDATA
            EdgeData *edge_data = edge_additions.E[i].edgeData;
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunction(source, destination, *edge_data,
                               vertex_value_old_next[source],
                               contrib_change_old, global_info_old);
          }
          if (ret) {
            if (use_lock) {
              vertex_locks[destination].writeLock();
              addToAggregation(contrib_change_old, delta[destination],
                               global_info_old);
              vertex_locks[destination].unlock();
            } else {
              addToAggregationAtomic(contrib_change_old, delta[destination],
                                     global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }

      parallel_for(long i = 0; i < edge_deletions.size; i++) {
        uintV source = edge_deletions.E[i].source;
        uintV destination = edge_deletions.E[i].destination;

        if (notDelZero(vertex_value_old_curr[source],
                      vertex_value_old_next[source], global_info_old) ||
            (forceActivateVertexForIteration(source, iter + 1,
                                             global_info_old))) {
          // Do retract for edge source->destination
          AggregationValueType contrib_change_old =
              aggregationValueIdentity<AggregationValueType>();

          // Do repropagate for edge source->destination.
          bool ret = false;
          if (use_delta_next_iteration) {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change_old, vertex_value_old_curr[source],
                vertex_value_old_next[source], global_info_old);
#ifdef EDGEDATA
            EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunctionDelta(source, destination, *edge_data,
                                    vertex_value_old_curr[source],
                                    vertex_value_old_next[source],
                                    contrib_change_old, global_info_old);
          } else {
            sourceChangeInContribution<AggregationValueType, VertexValueType,
                                       GlobalInfoType>(
                source, contrib_change_old,
                vertexValueIdentity<VertexValueType>(),
                vertex_value_old_next[source], global_info_old);
#ifdef EDGEDATA
            EdgeData *edge_data = edge_deletions.E[i].edgeData;
#else
            EdgeData *edge_data = &emptyEdgeData;
#endif
            ret = edgeFunction(source, destination, *edge_data,
                               vertex_value_old_next[source],
                               contrib_change_old, global_info_old);
          }
          if (ret) {
            if (use_lock) {
              vertex_locks[destination].writeLock();
              removeFromAggregation(contrib_change_old, delta[destination],
                                    global_info_old);
              vertex_locks[destination].unlock();

            } else {
              removeFromAggregationAtomic(contrib_change_old,
                                          delta[destination], global_info_old);
            }
            if (!changed[destination])
              changed[destination] = 1;
            if (!has_direct_changes)
              has_direct_changes = true;
          }
        }
      }
      phase_time = phase_timer.next();

      // Create frontier for next iteration
      vertexSubset temp_vs(n, frontier_curr);
      frontier_curr_vs = temp_vs;

      misc_time += phase_timer.next();
      iteration_time = iteration_timer.next();

      // Convergence check
      if (!has_direct_changes && frontier_curr_vs.isEmpty()) {
        if (iter == converged_iteration) {
          break;
        } else if (iter > converged_iteration) {
          assert(
              ("Missed switching to T-inc when iter == converged_iter", false));
        } else {
          // Values stable for the changed vertices at this iteration.
          // But, the changed vertices might receive new changes. So,
          // continue loop until iter == converged_iteration vertices may
          // still not have converged. So, keep continuing until
          // converged_iteration is reached.
        }
      }
      if (iter == 1) {
        iteration_time += pre_compute_time;
      }

      if (ae_enabled && shouldSwitch(iter, iteration_time)) {
        should_switch_now = true;
      }
      misc_time += phase_timer.stop();
      iteration_time += iteration_timer.stop();
    }

    cout << "Finished batch : " << full_timer.stop() << "\n";
    cout << "Number of iterations : " << converged_iteration << "\n";
    // testPrint();
    printOutput();
  }

  // Refactor this in a better way
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::my_graph;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::config;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::max_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::history_iterations;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::converged_iteration;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_lock;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_locks;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::aggregation_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_values;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::delta;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::use_source_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::source_change_in_contribution;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::n_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::global_info_old;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::vertex_value_old_prev;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::all;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_curr;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::frontier_next;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::changed;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::retract;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::propagate;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ingestor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::current_batch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::adaptive_executor;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::ae_enabled;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::testPrint;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::printOutput;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::shouldSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::performSwitch;
  using GraphBoltEngine<vertex, AggregationValueType, VertexValueType,
                        GlobalInfoType>::processVertexAddition;
};
#endif
