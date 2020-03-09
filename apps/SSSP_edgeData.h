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

#ifndef SSSP_EDGEDATA_H
#define SSSP_EDGEDATA_H

#include "../core/graph/edgeDataType.h"

/**
 *  Simple Edge Data Type for SSSP
 **/
struct SSSP_EdgeData : public EdgeDataType {
public:
  long weight = 0;
  SSSP_EdgeData() {}

  void createEdgeData(const char *edgeDataString) {
    weight = atol(edgeDataString);
  }

  void setEdgeDataFromPtr(EdgeDataType *edgeData) {
    weight = ((SSSP_EdgeData *)edgeData)->weight;
  }

  void del() {}

  std::string print() { return std::to_string(weight); }
};

// NOTE : The following typedef is important for the core files (graph.h, IO.h,
// etc.) to identify the type of the EdgeData.
// TODO : This is will be properly templatized in future.
typedef SSSP_EdgeData EdgeData;
#endif
