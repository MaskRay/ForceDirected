#ifndef FRUCHTERMANREINGOLD_HH
#define FRUCHTERMANREINGOLD_HH

#include <string.h>
#include <algorithm>
#include <array>
#include <functional>
#include <limits>
#include <vector>
#include <stdlib.h>
#include "Core.hh"
using std::function;
using std::min;
using std::max;
using std::vector;

template<typename T, size_t Dim>
struct FruchtermanReingold : ForceDirectedDrawing<T, Dim>
{
  FruchtermanReingold()
    : iterations(100) {}
  void set(const char* option, const char* value) {
    if (! strcmp(option, "iterations"))
      iterations = atoi(value);
  }
  virtual vector<Vector<T, Dim>> operator()(const Graph<T>& g) {
    vector<Vector<T, Dim>> pos(g.n), vel(g.n);
    T k = sqrt(accumulate(space.begin(), space.end(), T(1), std::multiplies<T>())) / g.n;
    function<T(T)> attractive = [&](T d) { return d * d / k; },
      repulsive = [&](T d) { return k * k / d; };
    for (int u = 0; u < g.n; u++) {
      T angle = 2 * M_PI * u / g.n;
      for (size_t dim = 0; dim < Dim; dim++)
        pos[u][dim] = space[dim] / 2;
      pos[u][0] += space[0] * cos(angle);
      pos[u][1] += space[1] * sin(angle);
    }
    for (int i = iterations; i > 0; i--) {
      T temperature = *std::min_element(space.begin(), space.end()) * i / iterations;
      for (int u = 0; u < g.n; u++)
        vel[u].fill(0);
      for (int u = 0; u < g.n; u++)
        for (int v = 0; v < g.n; v++)
          if (u != v) {
            Vector<T, Dim> dist = pos[v] - pos[u];
            vel[u] -= dist.unit() * repulsive(dist.norm());
          }
      for (int u = 0; u < g.n; u++)
        for (typename Graph<T>::It j = g.e[u].begin(); j != g.e[u].end(); j++)
          if (j->first != u) {
            Vector<T, Dim> dist = pos[j->first] - pos[u];
            vel[u] += dist.unit() * attractive(dist.norm());
          }
      for (int u = 0; u < g.n; u++) {
        pos[u] += vel[u].unit() * min(vel[u].norm(), temperature);
        for (size_t dim = 0; dim < Dim; dim++)
          pos[u][dim] = min(space[dim], max(T(0), pos[u][dim]));
      }
    }
    for (size_t dim = 0; dim < Dim; dim++) {
      T minx = std::numeric_limits<T>::max(),
        maxx = std::numeric_limits<T>::min();
      for (int u = 0; u < g.n; u++) {
        minx = min(minx, pos[u][dim]);
        maxx = max(maxx, pos[u][dim]);
      }
      for (int u = 0; u < g.n; u++)
        pos[u][dim] = (pos[u][dim] - minx) / (maxx - minx) * space[dim];
    }
    return pos;
  }
  array<T, Dim> space;
  int iterations;
};

#endif /* end of include guard: FRUCHTERMANREINGOLD_HH */
