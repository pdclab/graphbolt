// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
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

#ifndef UTIL_H
#define UTIL_H
#include "gettime.h"
#include "parseCommandLine.h"
#include "ligraUtils.h"
#include "parallel.h"
#include "rwlock.h"
#include <vector>
using namespace std;

#define ENABLE_CHECK_DESTINATION
// #define ENABLE_CHECK_SOURCE
#define CHECK_DEBUG_STATE(x) 1
#define ENABLE_PARTITION 1
#ifdef INCLUDEEXTRABUFFERSPACE
#define EDGE_ARRAY_BUFFER 10
#endif
#define USE_PARALLEL_SCAN

#define TIME_PRECISION 7
#define VAL_PRECISION 7
#define VAL_PRECISION2 13
#define PRINT_WIDTH 12

#define MOD_VAL 0.00190000001000000000

vector<int> debug_vertices = {624};

#endif // _UTIL_H
