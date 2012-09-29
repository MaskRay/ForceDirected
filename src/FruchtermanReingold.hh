#ifndef FRUCHTERMANREINGOLD_HH
#define FRUCHTERMANREINGOLD_HH

#include <string.h>
#include <vector>
#include <stdlib.h>
#include "Core.hh"
using std::vector;

class FruchtermanReingold : public ForceDirectedAlgorithm
{
public:
  FruchtermanReingold()
    : iterations(100) {}
  void set(const char* option, const char* value) {
    if (! strcmp(option, "iterations"))
      iterations = atoi(value);
  }
  virtual vector<Vector> operator()(const Graph&);
protected:
  int iterations;
};

#endif /* end of include guard: FRUCHTERMANREINGOLD_HH */
