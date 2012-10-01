#ifndef KAMADAKAWAI_HH
#define KAMADAKAWAI_HH

#include <string.h>
#include <stdlib.h>
#include "Core.hh"

template <size_t Size>
struct LinearSolver {};

template <>
struct LinearSolver<1>
{
  template<typename Vec>
    static Vec solve(double mat[1][1], Vec rhs) {
      return rhs / mat[0][0];
    }
};

// Cramer's rule
template <>
struct LinearSolver<2>
{
  template <typename Vec>
    static Vec solve(double mat[2][2], Vec rhs) {
      double denom = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
      double x_num = rhs[0]    * mat[1][1] - rhs[1]    * mat[0][1];
      double y_num = mat[0][0] * rhs[1]    - mat[1][0] * rhs[0]   ;
      Vec result;
      result[0] = x_num / denom;
      result[1] = y_num / denom;
      return result;
    }
};

template <>
struct LinearSolver<3>
{
  template <typename Vec>
    static Vec solve(double mat[2][2], Vec rhs) {
      double denom = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
        - mat[1][0] * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2])
        + mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
      double x_num = rhs[0]    * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
        - rhs[1]    * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2])
        + rhs[2]    * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
      double y_num = mat[0][0] * (rhs[1]    * mat[2][2] - rhs[2]    * mat[1][2])
        - mat[1][0] * (rhs[0]    * mat[2][2] - rhs[2]    * mat[0][2])
        + mat[2][0] * (rhs[0]    * mat[1][2] - rhs[1]    * mat[0][2]);
      double z_num = mat[0][0] * (mat[1][1] * rhs[2]    - mat[2][1] * rhs[1]   )
        - mat[1][0] * (mat[0][1] * rhs[2]    - mat[2][1] * rhs[0]   )
        + mat[2][0] * (mat[0][1] * rhs[1]    - mat[1][1] * rhs[0]   );
      Vec result;
      result[0] = x_num / denom;
      result[1] = y_num / denom;
      result[2] = z_num / denom;
      return result;
    }
};

template<typename T, size_t Dim>
struct KamadaKawai : ForceDirectedDrawing<T, Dim>
{
  using ForceDirectedDrawing<T, Dim>::space;

  KamadaKawai(const array<T, Dim>& space)
    : ForceDirectedDrawing<T, Dim>(space)
      , tolerance(1e-6)
      , spring_strength(1) {}
  void set(const char* option, const char* value) {
    if (! strcmp(option, "tolerance"))
      tolerance = atof(value);
  }
  virtual void operator()(const Graph<T>& g, vector<Vector<T, Dim>>& pos) {
    vector<vector<T>> dist(g.n, vector<T>(g.n, numeric_limits<T>::max())),
      strength(g.n, vector<T>(g.n));
    for (int u = 0; u < g.n; u++)
      for (typename Graph<T>::It j = g.e[u].begin(); j != g.e[u].end(); j++)
        dist[u][j->first] = dist[j->first][u] = j->second;
    for (int u = 0; u < g.n; u++)
      dist[u][u] = T(0);
    // Floyd-Warshall
    // TODO resort to Johnson's in a sparse graph
    for (int p = 0; p < g.n; p++)
      for (int u = 0; u < g.n; u++)
        for (int v = 0; v < g.n; v++)
          dist[u][v] = min(dist[u][v], dist[u][p] + dist[p][v]);

    double edge_length = 0;
    for (int u = 0; u < g.n; u++)
      for (int v = u; ++v < g.n; )
        edge_length = max(edge_length, dist[u][v]);
    edge_length = *std::max_element(space.begin(), space.end()) / edge_length;
    for (int u = 0; u < g.n; u++)
      for (int v = u; ++v < g.n; ) {
        T d = dist[u][v];
        dist[u][v] = dist[v][u] = edge_length * d;
        strength[u][v] = strength[v][u] = spring_strength / (d * d);
      }

    // contribution of vertex v to vertex u
    function<Vector<T, Dim>(int, int)> compute_partial_deriv = [&](int u, int v) {
      if (u == v) {
        Vector<T, Dim> res = {};
        res.fill(T(0));
        return res;
      }
      Vector<T, Dim> diff = pos[u] - pos[v];
      T l = dist[u][v];
      T d = diff.norm();
      return (diff - diff * (l / d)) * strength[u][v];
    };
    function<Vector<T, Dim>(int)> compute_partial_derivs = [&](int u) {
      Vector<T, Dim> res = {};
      for (int v = 0; v < g.n; v++)
        res += compute_partial_deriv(u, v);
      return res;
    };

    // find max delta
    LayoutTolerance done(tolerance);
    int pivot = 0;
    T max_delta(0);
    vector<Vector<T, Dim>> partials(g.n);
    for (int u = 0; u < g.n; u++) {
      partials[u] = compute_partial_derivs(u);
      T delta = partials[u].norm();
      if (delta > max_delta) {
        pivot = u;
        max_delta = delta;
      }
    }

    while (! done(max_delta, true)) {
      vector<Vector<T, Dim>> p_partials(g.n);
      for (int u = 0; u < g.n; u++)
        p_partials[u] = compute_partial_deriv(u, pivot);
      // tune vertex pivot
      do {
        //double E = 0;
        //for (int u = 0; u < g.n; u++)
          //for (int v = u; ++v < g.n; ) {
            //double l = dist[u][v],
                   //d = (pos[u] - pos[v]).norm();
            //E += 0.5 * strength[u][v] * pow(d-l, 2);
          //}
        //fprintf(stderr, "E = %lf\n", E);

        T ddE[Dim][Dim] = {};
        for (int u = 0; u < g.n; u++)
          if (pivot != u) {
            Vector<T, Dim> diff = pos[pivot] - pos[u];
            T d2 = diff.norm2(), d = sqrt(d2), inv_d3 = T(1) / (d2 * d),
              l = dist[pivot][u],
              k = strength[pivot][u];
            for (size_t i = 0; i < Dim; i++)
              for (size_t j = 0; j < Dim; j++)
                if (i == j)
                  ddE[i][j] += k * (T(1) + (l * (diff[i] * diff[i] - d2) * inv_d3));
                else
                  ddE[i][j] += k * l * diff[i] * diff[j] * inv_d3;
          }

        Vector<T, Dim> step = LinearSolver<Dim>::solve(ddE, - partials[pivot]);
        for (size_t dim = 0; dim < Dim; dim++)
          pos[pivot][dim] += step[dim];
        partials[pivot] = compute_partial_derivs(pivot);
        max_delta = partials[pivot].norm();
      } while (! done(max_delta, false));

      int old_p = pivot;
      for (int u = 0; u < g.n; u++) {
        Vector<T, Dim> old = p_partials[u],
          new_ = compute_partial_deriv(u, old_p);
        partials[u] += new_ - old;
        T delta = partials[u].norm();
        if (delta > max_delta) {
          pivot = u;
          max_delta = delta;
        }
      }
    }

    normalizeToSpace(pos, space);
  }

  struct LayoutTolerance
  {
    LayoutTolerance(T tolerance)
      : tolerance(tolerance)
      , last(numeric_limits<T>::max())
      , last_l(numeric_limits<T>::max()) {}
    bool operator()(T delta, bool global) {
      if (global) {
        if (last == numeric_limits<T>::max()) {
          last = delta;
          return false;
        }
        T diff = last - delta;
        if (diff < 0) diff = - diff;
        last = delta;
        return delta < tolerance;
      } else {
        if (last_l == numeric_limits<T>::max()) {
          last_l = delta;
          return false;
        }
        T diff = last_l - delta;
        if (diff < 0) diff = - diff;
        last_l = delta;
        return delta < tolerance;
      }
    }
    T tolerance, last, last_l;
  };

  T tolerance, spring_strength;
};

#endif /* end of include guard: KAMADAKAWAI_HH */
