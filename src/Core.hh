#ifndef CORE_HH
#define CORE_HH

#include <array>
#include <cmath>
#include <utility>
#include <vector>
using std::array;
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
  virtual ~ForceDirectedDrawing() {}
  virtual std::vector<Vector<T, Dim> > operator()(const Graph<T>&) = 0;
};

#endif /* end of include guard: CORE_HH */
