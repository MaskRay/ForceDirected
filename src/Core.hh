#ifndef CORE_HH
#define CORE_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <functional>
#include <utility>
#include <vector>
using std::array;
using std::bind;
using std::function;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::sqrt;
using std::vector;

template<typename T, size_t Dim>
struct Vector : array<T, Dim>
{
  T dot(const Vector& r) const {
    T s = 0;
    for (int i = 0; i < Dim; i++) s += (*this)[i] * r[i];
    return s;
  }
  T dist(const Vector& r) const { return (*this - r).norm(); }
  Vector unit() const { return (*this) / norm(); }
  Vector operator-() const {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; i++) res[i] = - (*this)[i];
    return res;
  }
  Vector operator*(T r) const {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; i++) res[i] = (*this)[i] * r;
    return res;
  }
  Vector operator/(T r) const {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; i++) res[i] = (*this)[i] / r;
    return res;
  }
  Vector operator+(const Vector& r) const {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; i++) res[i] = (*this)[i] + r[i];
    return res;
  }
  Vector operator-(const Vector& r) const {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; i++) res[i] = (*this)[i] - r[i];
    return res;
  }
  Vector operator+=(const Vector& r) {
    for (int i = 0; i < Dim; i++) (*this)[i] += r[i];
    return *this;
  }
  Vector operator-=(const Vector& r) {
    for (int i = 0; i < Dim; i++) (*this)[i] -= r[i];
    return *this;
  }
  T norm() const { return sqrt(norm2()); }
  T norm2() const {
    T s(0);
    for (int i = 0; i < Dim; i++) s += (*this)[i] * (*this)[i];
    return s;
  }
};

template<typename T, size_t Dim>
Vector<T, Dim> operator*(T l, const Vector<T, Dim>& r)
{ return r * l; }

template<typename T>
class Graph
{
public:
  Graph(int n) : n(n), e(n) {}
  void addEdge(int u, int v, const T& w) {
    e[u].push_back(std::make_pair(v, w));
    e[v].push_back(std::make_pair(u, w));
  }
  typedef typename vector<pair<int, T> >::const_iterator It;
  int n;
  vector<vector<pair<int, T>>> e;
};

template<typename T, size_t Dim>
struct ForceDirectedDrawing
{
  ForceDirectedDrawing(const array<T, Dim>& space) : space(space) {}
  virtual ~ForceDirectedDrawing() {}
  virtual void operator()(const Graph<T>&, vector<Vector<T, Dim>>&) = 0;
  array<T, Dim> space;
};

template<typename T, size_t Dim, typename G>
void normalizeToSpace(vector<Vector<T, Dim>>& pos, const G& space)
{
  for (size_t dim = 0; dim < Dim; dim++) {
    T minx = numeric_limits<T>::max(),
      maxx = numeric_limits<T>::min();
    for (auto &i : pos) {
      minx = min(minx, i[dim]);
      maxx = max(maxx, i[dim]);
    }
    for (auto &i : pos)
      i[dim] = (i[dim] - minx) / (maxx - minx) * space[dim];
  }
}

#endif /* end of include guard: CORE_HH */
