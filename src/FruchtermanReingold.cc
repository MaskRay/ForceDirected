#include "FruchtermanReingold.hh"
#include <algorithm>
#include <functional>
using namespace std;

vector<Vector> FruchtermanReingold::operator()(const Graph& g)
{
  vector<Vector> pos(g.n), vel(g.n);
  ft k = sqrt(width * height / g.n);
  function<ft(ft)> attractive = [&](ft d) { return d * d / k; },
    repulsive = [&](ft d) { return k * k / d; };
  for (int u = 0; u < g.n; u++) {
    pos[u].x = width / 2 + min(width, height) * cos(2 * M_PI * u / g.n);
    pos[u].y = height / 2 + min(width, height) * sin(2 * M_PI * u / g.n);
  }
  for (int i = 10; i > 0; i--) {
    ft temperature = min(width, height) * i / iterations;
    std::fill(vel.begin(), vel.end(), Vector(0.0, 0.0));
    for (int u = 0; u < g.n; u++)
      for (int v = 0; v < g.n; v++)
        if (u != v) {
          Vector dist = pos[v] - pos[u];
          vel[u] -= dist.unit() * repulsive(dist.norm());
        }
    for (int u = 0; u < g.n; u++)
      for (Graph::It j = g.e[u].begin(); j != g.e[u].end(); j++)
        if (*j != u) {
          Vector dist = pos[*j] - pos[u];
          vel[u] += dist.unit() * attractive(dist.norm());
        }
    for (int u = 0; u < g.n; u++) {
      pos[u] += vel[u].unit() * min(vel[u].norm(), temperature);
      pos[u].x = min(width, max(0.0, pos[u].x));
      pos[u].y = min(height, max(0.0, pos[u].y));
    }
  }
  return pos;
}
