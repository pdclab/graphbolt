#ifndef CF_EDGEDATA_H
#define CF_EDGEDATA_H

#include "../core/graph/edgeDataType.h"

/**
 *  Simple Edge Data Type for CF
 **/
struct CF_EdgeData : public EdgeDataType {
public:
  double weight = 0;
  CF_EdgeData() {}

  void createEdgeData(const char *edgeDataString) {
    weight = atof(edgeDataString);
  }

  void setEdgeDataFromPtr(EdgeDataType *edgeData) {
    weight = ((CF_EdgeData *)edgeData)->weight;
  }

  void del() {}

  std::string print() { return std::to_string(weight); }
};

typedef CF_EdgeData EdgeData;
#endif