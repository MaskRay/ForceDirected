#ifndef KAMADAKAWAI_HH
#define KAMADAKAWAI_HH

#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <vector>
#include "Core.hh"

template<typename T>
class KamadaKawai : public ForceDirectedAlgorithm<T>
{
  using std::bind;
  using std::vector;
public:
  KamadaKawai()
    : tolerance(100) {}
  void set(const char* option, const char* value) {
    if (! strcmp(option, "tolerance"))
      tolerance = atoi(value);
  }
  virtual vector<Vector> operator()(const Graph& g) {
    vector<vector<T>> dist(n, vector(n, numeric_limits<T>::max()));
    // Floyd-Warshall
    // TODO resort to Johnson's in a sparse graph
    for (int p = 0; p < g.n; p++)
      for (int u = 0; u < g.n; u++)
        for (int v = 0; v < g.n; v++)
          dist[u][v] = min(dist[u][v], dist[u][p] + dist[p][v]);

    // spring strength
    function<T(int)> compute_strength = [&](T dist) {
      return spring_strength / (dist * dist);
    }

    // contribution of vertex v to vertex u
    function<T(int)> compute_partial_deriv = [&](int u, int v) {
      if (u == v)
        return Vector<T>(0, 0);
      Vector d = pos[u] - pos[v];
      return compute_strength(dist[u][v]) * (d - d.unit() * dist[u][v]);
    };
    function<T(int)> compute_partial_derivs = [&](int u) {
      for (int v = 0; v < g.n; v++)
        res += compute_partial_deriv(u, v);
      return res;
    };

    // find max delta
    LayoutTolerance done(tolerance);
    int pivot = 0;
    T max_delta = 0;
    vector<Vector<T>> partial_derivs(g.n);
    for (int u = 0; u < g.n; u++) {
      Vector<T> deriv = compute_partial_derivs(u);
      partial_derivs[u] = deriv;
      T delta = deriv.norm();
      if (delta > max_delta) {
        pivot = u;
        max_delta = delta;
      }
    }

    while (! done(max_delta, true)) {
      vector<T> partials(g.n);
      for (int u = 0; u < g.n; u++)
        partials[u] = compute_partial_deriv(pivot, u);
      // tune vertex pivot
      do {
        double E = 0;
        for (int u = 0; u < g.n; u++)
          for (int v = 0; v < g.n; v++) {
            double d = dist[u][v],
                   k = compute_strength(d);
            E += 0.5 * k * pow(d-k, 2);
          }
        printf("E = %lf\n", E);

        vector<vector<T>> ddE(2, 2);
      } while (! done(delta, false));
    }
  }
protected:
  T tolerance, spring_strength;
  struct LayoutTolerance
  {
    LayoutTolerance(T tolerance = T(1e-3))
      : tolerance(tolerance)
      , last(numeric_limits<T>::max())
      , last_l(numeric_limits<T>::max())
    bool operator()(T delta, bool global) {
      if (global) {
        if (last == std::numeric_limits<T>::max()) {
          last = delta;
          return false;
        }
        T diff = last - delta;
        if (diff < 0) diff = - diff;
        bool done = diff == T(0) || diff / last < tolerance;
        last = delta;
        return done;
      } else {
        if (last_l == std::numeric_limits<T>::max()) {
          last_l = delta;
          return false;
        }
        T diff = last_l - delta;
        if (diff < 0) diff = - diff;
        bool done = diff == T(0) || diff / last_l < tolerance;
        last_l = delta;
        return done;
      }
    }
    T tolerance;
  };
};

#endif /* end of include guard: KAMADAKAWAI_HH */
