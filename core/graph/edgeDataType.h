#ifndef EDGEDATATYPE_H
#define EDGEDATATYPE_H

#include "../common/utils.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

struct EdgeDataType {
public:
  virtual void createEdgeData(const char *edgeDataString) = 0;
  virtual void setEdgeDataFromPtr(EdgeDataType *edgeData) = 0;
  // virtual void print(std::ostream& os) = 0;
  virtual std::string print() = 0;
  virtual void del() = 0;
  virtual ~EdgeDataType() = default;
};

#endif