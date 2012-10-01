#ifndef WALSHAW_HH
#define WALSHAW_HH

#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "Core.hh"
#include "KdTree.hh"

template<typename T, size_t Dim>
struct Walshaw : ForceDirectedDrawing<T, Dim>
{
  using ForceDirectedDrawing<T, Dim>::space;

  Walshaw(const array<T, Dim>& space)
    : ForceDirectedDrawing<T, Dim>(space)
      , use_BSP(true)
      , separation_constant(2)
      , force_constant(0.01)
      , alpha(1.5)
      , iterations(50) {}
  virtual void operator()(const Graph<T>& g, vector<Vector<T, Dim>>& pos) {
    vector<Vector<T, Dim>> vel(g.n);
    T k = separation_constant * std::pow(accumulate(space.begin(), space.end(), T(1), std::multiplies<T>()) / T(g.n), T(1) / T(Dim));
    function<T(T)> global = [&](T d) { return - k * k / d * force_constant; };

    for (int i = iterations; i > 0; i--) {
      T temperature = *std::min_element(space.begin(), space.end()) * i / iterations;
      for (int u = 0; u < g.n; u++)
        vel[u].fill(0);
      if (use_BSP) {
        KdTree<T, Dim> kd(pos, alpha);
        for (int u = 0; u < g.n; u++)
          vel[u] += kd.getRepulsive(pos[u]) * k * k * force_constant;
      } else {
        for (int u = 0; u < g.n; u++)
          for (int v = 0; v < g.n; v++)
            if (u != v) {
              Vector<T, Dim> dist = pos[v] - pos[u];
              vel[u] += dist.unit() * global(dist.norm());
            }
      }
      for (int u = 0; u < g.n; u++)
        for (typename Graph<T>::It j = g.e[u].begin(); j != g.e[u].end(); j++)
          if (j->first != u) {
            Vector<T, Dim> dist = pos[j->first] - pos[u];
            vel[u] += dist.unit() * ((dist.norm() - k) / g.e[u].size() - global(dist.norm()));
          }
      for (int u = 0; u < g.n; u++)
        pos[u] += vel[u].unit() * min(vel[u].norm(), temperature);
    }

    normalizeToSpace(pos, space);
  }

  int iterations;
  T separation_constant, force_constant;

  // if BSP
  bool use_BSP;
  T alpha;
};

#endif /* end of include guard: WALSHAW_HH */
