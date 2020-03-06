#ifndef LP_EDGEDATA_H
#define LP_EDGEDATA_H

#include "../core/graph/edgeDataType.h"

/**
 *  Complex Edge Data Type for LP
 **/
struct LP_EdgeData : public EdgeDataType {
public:
  double weight = 0;
  LP_EdgeData() {}

  void createEdgeData(const char *edgeDataString) {
    weight = atof(edgeDataString);
  }

  void setEdgeDataFromPtr(EdgeDataType *edgeData) {
    weight = ((LP_EdgeData *)edgeData)->weight;
  }

  void del() {}

  std::string print() { return std::to_string(weight); }
};

typedef LP_EdgeData EdgeData;

#endif
