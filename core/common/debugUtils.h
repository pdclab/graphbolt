#ifndef DEBUG_UTIL_H
#define DEBUG_UTIL_H

#include <vector>
#include "parallel.h"
using namespace std;

#define ENABLE_TIMER_LOGS 1
#ifndef ENABLE_TIMER_LOGS
#define DEBUG_LOG_1 1
#define DEBUG_LOG_2 1
// #define DEBUG_LOG_3 1
// #define DEBUG_LOG_4 1
// #define DEBUG_LOG_5 1
#define OUTPUT_LOG 1
#define PRINT_HISTORY 1
#endif
// #define MOD_VAL 0.20000001
#define DEBUG_FEATURE 1
#define DEBUG_STATE 1

#define ENABLE_CHECK_DESTINATION

vector<uintV> debug_vertices = {2481};
vector<uintV> debugVertex1 = {570};

#ifdef ENABLE_CHECK_DESTINATION
#define CHECK_DESTINATION(x)                                                   \
  std::find(debug_vertices.begin(), debug_vertices.end(), x) != debug_vertices.end()
#else
#define CHECK_DESTINATION(x) false
// #define CHECK_DESTINATION(x) true
#endif

#ifdef ENABLE_CHECK_SOURCE
#define CHECK_SOURCE(x)                                                        \
  std::find(debugVertex1.begin(), debugVertex1.end(), x) != debugVertex1.end()
#else
#define CHECK_SOURCE(x) true
#endif

#define DEBUG_MODE 0
#define DEBUG_LEVEL 11
#define VERTEX_HISTORY 8
#define VERTEX_DEBUG_LEVEL 9
#define EDGE_DEBUG_LEVEL 10
#define EDGE_DETAILED_DEBUG_LEVEL 11

#if DEBUG_MODE
#define MY_DEBUG_LOG(level, body)                                              \
  {                                                                            \
    if (level <= DEBUG_LEVEL) {                                                \
      body();                                                                  \
    }                                                                          \
  }
#else
#define MY_DEBUG_LOG(level, body)                                              \
  {}
#endif

#if DEBUG_MODE
#else
#define TIMER_LOGS 1
#endif
#if TIMER_LOGS
#define MY_TIMER_LOGS(body)                                                    \
  { body(); }
#else
#define MY_TIMER_LOGS(body)                                                    \
  {}
#endif

#endif
