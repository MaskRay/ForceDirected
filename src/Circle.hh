#ifndef CIRCLE_HH
#define CIRCLE_HH

#include <cmath>
#include "Core.hh"

template<typename T, size_t Dim>
struct Circle : ForceDirectedDrawing<T, Dim>
{
  using ForceDirectedDrawing<T, Dim>::space;

  Circle(const array<T, Dim>& space)
    : ForceDirectedDrawing<T, Dim>(space) {}
  virtual void operator()(const Graph<T>& g, vector<Vector<T, Dim>>& pos) {
    for (int u = 0; u < g.n; u++) {
      T angle = 2 * M_PI * u / g.n;
      for (size_t dim = 0; dim < Dim; dim++)
        pos[u][dim] = space[dim] / 2;
      pos[u][0] += space[0] * std::cos(angle);
      pos[u][1] += space[1] * std::sin(angle);
    }
  }
};

#endif /* end of include guard: CIRCLE_HH */
