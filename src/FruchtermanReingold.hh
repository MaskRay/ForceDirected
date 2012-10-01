#ifndef FRUCHTERMANREINGOLD_HH
#define FRUCHTERMANREINGOLD_HH

#include <string.h>
#include <stdlib.h>
#include "Core.hh"
#include "KdTree.hh"

template<typename T, size_t Dim>
struct FruchtermanReingold : ForceDirectedDrawing<T, Dim>
{
  using ForceDirectedDrawing<T, Dim>::space;

  FruchtermanReingold(const array<T, Dim>& space)
    : ForceDirectedDrawing<T, Dim>(space)
      , use_BSP(true)
      , alpha(1.5)
      , iterations(50) {}
  virtual void operator()(const Graph<T>& g, vector<Vector<T, Dim>>& pos) {
    vector<Vector<T, Dim>> vel(g.n);
    T k = sqrt(accumulate(space.begin(), space.end(), T(1), std::multiplies<T>())) / g.n;
    function<T(T)> attractive = [&](T d) { return d * d / k; },
      repulsive = [&](T d) { return k * k / d; };

    for (int i = iterations; i > 0; i--) {
      T temperature = *std::min_element(space.begin(), space.end()) * i / iterations;
      for (int u = 0; u < g.n; u++)
        vel[u].fill(0);
      if (use_BSP) {
        KdTree<T, Dim> kd(pos, alpha);
        for (int u = 0; u < g.n; u++)
          vel[u] += kd.getRepulsive(pos[u]) * k * k;
      } else {
        for (int u = 0; u < g.n; u++)
          for (int v = 0; v < g.n; v++)
            if (u != v) {
              Vector<T, Dim> dist = pos[v] - pos[u];
              vel[u] -= dist.unit() * repulsive(dist.norm());
            }
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

    normalizeToSpace(pos, space);
  }

  int iterations;

  // if BSP
  bool use_BSP;
  T alpha;
};

#endif /* end of include guard: FRUCHTERMANREINGOLD_HH */
