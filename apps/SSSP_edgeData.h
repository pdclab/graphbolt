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

typedef SSSP_EdgeData EdgeData;
#endif