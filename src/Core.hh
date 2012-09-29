#ifndef CORE_HH
#define CORE_HH

#include <cmath>
#include <vector>
using std::sqrt;

#if defined(USE_FLOAT)
typedef float ft;
#else
typedef double ft;
#endif

struct Vector
{
  ft x, y;
  Vector() {}
  Vector(ft x, ft y) : x(x), y(y) {}
  ft dot(const Vector& r) const {
    return x * r.x + y * r.y;
  }
  ft dist(const Vector& r) const {
    return (*this - r).norm();
  }
  Vector unit() const {
    ft n = 1.0 / norm();
    return Vector(x * n, y * n);
  }
  Vector operator*(ft r) const { return Vector(x * r, y * r); }
  Vector operator+(const Vector& r) const { return Vector(x + r.x, y + r.y); }
  Vector operator-(const Vector& r) const { return Vector(x - r.x, y - r.y); }
  Vector operator+=(const Vector& r) {
    x += r.x;
    y += r.y;
    return *this;
  }
  Vector operator-=(const Vector& r) {
    x -= r.x;
    y -= r.y;
    return *this;
  }
  ft norm() const { return sqrt(x * x + y * y); }
  ft norm2() const { return x * x + y * y; }
};

Vector operator*(ft l, const Vector& r);

class Graph
{
public:
  Graph(int n) : n(n), e(n) {}
  void addEdge(int u, int v) {
    e[u].push_back(v);
    e[v].push_back(u);
  }
  typedef std::vector<int>::const_iterator It;
  int n;
  std::vector<std::vector<int> > e;
};

class ForceDirectedAlgorithm
{
public:
  virtual ~ForceDirectedAlgorithm() {}
  virtual std::vector<Vector> operator()(const Graph&) = 0;
  void setSize(ft width, ft height) { this->width = width; this->height = height; }
protected:
  ft width, height;
};

ft getRand(ft);

#endif /* end of include guard: CORE_HH */
