#ifndef KDTREE_HH
#define KDTREE_HH

#include <vector>
#include "Core.hh"

template<typename T, size_t Dim>
class KdTree
{
public:
  KdTree(const vector<Vector<T, Dim>>& coords, T alpha)
    : coords(coords)
      , alpha(alpha)
      , node_size_threashold(4) {
      root = build(0, coords.size(), 0);
    }
  ~KdTree() { delete root; }
  struct Node {
    Cube<T, Dim> bounding; // minimum bounding box
    Vector<T, Dim> sum; // sum of coordinates in subtree
    int L, R; // coords[L..R)
    Node* ch[2];
    Node() { ch[0] = ch[1] = NULL; }
    ~Node() { delete ch[0]; delete ch[1]; }
    bool isLeaf() { return ch[0] == NULL && ch[1] == NULL; }
  };
  Vector<T, Dim> getRepulsive(const Vector<T, Dim>& orig) {
    return getRepulsive(root, orig);
  }
protected:
  Node* build(int L, int R, size_t dim) {
    if (L == R) return NULL;
    Node* rt = new Node;
    rt->L = L;
    rt->R = R;
    if (R - L <= node_size_threashold) {
      rt->sum.fill(0);
      rt->bounding.lo = rt->bounding.hi = coords[L];
      for (int i = L; i < R; i++) {
        rt->sum += coords[i];
        for (size_t dim = 0; dim < Dim; dim++) {
          rt->bounding.lo[dim] = min(rt->bounding.lo[dim], coords[i][dim]);
          rt->bounding.hi[dim] = max(rt->bounding.hi[dim], coords[i][dim]);
        }
      }
    } else {
      int M = (L + R) / 2;
      std::nth_element(coords.begin() + L, coords.begin() + M, coords.begin() + R,
          [&](const Vector<T, Dim>& l, const Vector<T, Dim>& r) { return l[dim] < r[dim]; });
      T median = coords[M][dim];

      rt->ch[0] = build(L, M, (dim + 1) % Dim);
      rt->ch[1] = build(M, R, (dim + 1) % Dim);

      rt->sum = rt->ch[0]->sum + rt->ch[1]->sum;
      rt->bounding = rt->ch[0]->bounding;
      for (size_t dim = 0; dim < Dim; dim++) {
        rt->bounding.lo[dim] = min(rt->bounding.lo[dim], rt->ch[1]->bounding.lo[dim]);
        rt->bounding.hi[dim] = max(rt->bounding.hi[dim], rt->ch[1]->bounding.hi[dim]);
      }
    }
    return rt;
  }
  Vector<T, Dim> getRepulsive(Node* rt, const Vector<T, Dim>& orig) {
    if (rt == NULL) {
      Vector<T, Dim> res;
      res.fill(0);
      return res;
    }
    if (rt->isLeaf()) {
      Vector<T, Dim> res;
      res.fill(0);
      for (int i = rt->L; i < rt->R; i++)
        if (orig != coords[i]) { // exclude itself
          Vector<T, Dim> diff = orig - coords[i];
          res += diff.unit() / diff.norm();
        }
      return res;
    }

    // approximate
    T measure(0);
    Vector<T, Dim> barycenter = rt->sum / (rt->R - rt->L);
    for (size_t dim = 0; dim < Dim; dim++)
      measure = max(measure, rt->bounding.hi[dim] - rt->bounding.lo[dim]);
    auto diff = orig - barycenter;
    T d = diff.norm();
    if (d / measure > alpha)
      return diff.unit() / d * (rt->R - rt->L);

    return getRepulsive(rt->ch[0], orig) + getRepulsive(rt->ch[1], orig);
  }

  T alpha;
  Node* root;
  vector<Vector<T, Dim>> coords;
  size_t node_size_threashold;
};

#endif /* end of include guard:  */
