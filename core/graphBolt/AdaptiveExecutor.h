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

#ifndef __ADAPTIVE_EXECUTOR_H__
#define __ADAPTIVE_EXECUTOR_H__

#include "../common/utils.h"
#include "../graph/graph.h"
#include "../graph/vertexSubset.h"
#include <limits>

class AdaptiveExecutor {
  int history_iterations;
  int converged_iteration;

  long *active_edges;
  double *edge_map_time;
  double *vertex_map_time;
  double *copy_time;
  double *misc_time;

  double avg_active_edges;
  double avg_edge_map_time;
  double avg_vertex_map_time;
  double avg_copy_time;
  double avg_misc_time;

  double edge_map_slope;
  double edge_map_intercept;

  double switch_cost;
  double approximate_time_for_curr_iter;

public:
  // Iterations range from 1 to history_iterations. So, set
  // t_history_iterations+1 as the history_iterations
  AdaptiveExecutor(int t_history_iterations)
      : history_iterations(t_history_iterations + 1), converged_iteration(100) {
    if (history_iterations > 0) {
      active_edges = new long[history_iterations];
      edge_map_time = new double[history_iterations];
      copy_time = new double[history_iterations];
      vertex_map_time = new double[history_iterations];
      misc_time = new double[history_iterations];
    }
    reset();
  }

  void reset() {
    for (int i = 0; i < history_iterations; i++) {
      active_edges[i] = 0.0;
      edge_map_time[i] = 0.0;
      copy_time[i] = 0.0;
      vertex_map_time[i] = 0.0;
      misc_time[i] = 0.0;
    }

    avg_active_edges = 0.0;
    avg_edge_map_time = 0.0;
    avg_vertex_map_time = 0.0;
    avg_copy_time = 0.0;
    avg_misc_time = 0.0;

    edge_map_slope = 0.0;
    edge_map_intercept = 0.0;

    switch_cost = 0.0;
    approximate_time_for_curr_iter = numeric_limits<double>::max();
  }

  ~AdaptiveExecutor() {
    if (history_iterations != 0) {
      delete[] active_edges;
      delete[] edge_map_time;
      delete[] vertex_map_time;
      delete[] copy_time;
      delete[] misc_time;
    }
  }

  template <class vertex>
  inline void updateEdgesProcessed(int iter, graph<vertex> &t_graph,
                                   vertexSubset &frontier) {
    long edges_to_process =
        sequence::plusReduceDegree(t_graph.V, frontier.d, (long)t_graph.n);
    active_edges[iter] = edges_to_process;
  }

  inline void updateEdgeMapTime(int iter, double t_edge_map_time) {
    edge_map_time[iter] = t_edge_map_time;
  }

  inline void updateVertexMapTime(int iter, double t_vertex_map_time) {
    vertex_map_time[iter] = t_vertex_map_time;
  }

  inline void updateCopyTime(int iter, double t_copy_time) {
    copy_time[iter] = t_copy_time;
  }

  inline void updateMiscTime(int iter, double t_misc_time) {
    misc_time[iter] = t_misc_time;
  }

  inline double approximateTimeForCurrIter() {
    return approximate_time_for_curr_iter;
  }

  inline void setApproximateTimeForCurrIter(double t_curr_iter_time) {
    approximate_time_for_curr_iter = t_curr_iter_time;
  }

  inline void updateApproximateTimeForEdges(long t_active_edges) {
    approximate_time_for_curr_iter = t_active_edges * edge_map_slope +
                                     edge_map_intercept + avg_vertex_map_time +
                                     avg_copy_time + avg_misc_time;
  }

  // TODO : Rename function
  void updateEquation(int t_converged_iter) {
    converged_iteration = t_converged_iter;

    avg_active_edges = 0.0;
    avg_edge_map_time = 0.0;
    avg_vertex_map_time = 0.0;
    avg_copy_time = 0.0;
    avg_misc_time = 0.0;
    for (int i = 1; i <= converged_iteration; i++) {
      avg_active_edges += active_edges[i];
      avg_edge_map_time += edge_map_time[i];
      avg_vertex_map_time += vertex_map_time[i];
      avg_copy_time += copy_time[i];
      avg_misc_time += misc_time[i];
    }
    avg_active_edges = (double)avg_active_edges / (double)converged_iteration;
    avg_edge_map_time = avg_edge_map_time / converged_iteration;
    avg_vertex_map_time = avg_vertex_map_time / converged_iteration;
    avg_copy_time = avg_copy_time / converged_iteration;
    avg_misc_time = avg_misc_time / converged_iteration;

    double sum_top = 0.0;
    double sum_bottom = 0.0;
    double temp = 0.0;
    for (int i = 1; i <= converged_iteration; i++) {
      temp = active_edges[i] - avg_active_edges;
      sum_top += temp * (edge_map_time[i] - avg_edge_map_time);
      sum_bottom += (temp * temp);
    }
    edge_map_slope = sum_top / sum_bottom;
    edge_map_intercept = avg_edge_map_time - edge_map_slope * avg_active_edges;
  }

  void print() {
    cout << "Averages, " << setprecision(TIME_PRECISION) << edge_map_slope
         << ", " << edge_map_intercept << ", " << avg_vertex_map_time << ", "
         << avg_copy_time << ", " << avg_misc_time << endl;
  }
};

#endif
