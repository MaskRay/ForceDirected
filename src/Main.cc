#include "Core.hh"
#include "Circle.hh"
#include "FruchtermanReingold.hh"
#include "Walshaw.hh"
#include "KamadaKawai.hh"
#include "Cmdline.h"

int main(int argc, char* argv[])
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0)
    return 1;

  bool use_w = false;
  ForceDirectedDrawing<double, 2>* algo = NULL;
  array<double, 2> space;
  switch (args_info.algorithm_arg) {
  case 0:
    {
      space[0] = args_info.x_arg;
      space[1] = args_info.y_arg;
      auto a = new FruchtermanReingold<double, 2>(space);
      a->iterations = args_info.iterations_arg;
      a->separation_constant = args_info.separation_arg;
      a->force_constant = args_info.repulsive_arg;
      if (a->use_BSP = args_info.kd_arg != 0)
        a->alpha = args_info.alpha_arg;
      algo = a;
    }
    break;
  case 1:
    {
      space[0] = args_info.x_arg;
      space[1] = args_info.y_arg;
      auto a = new Walshaw<double, 2>(space);
      a->iterations = args_info.iterations_arg;
      a->separation_constant = args_info.separation_arg;
      a->force_constant = args_info.repulsive_arg;
      if (a->use_BSP = args_info.kd_arg != 0)
        a->alpha = args_info.alpha_arg;
      algo = a;
    }
    break;
  case 2:
    {
      space[0] = args_info.x_arg;
      space[1] = args_info.y_arg;
      auto a = new KamadaKawai<double, 2>(space);
      algo = a;
      use_w = true;
    }
    break;
  default:
    return 2;
  }

  int n, m;
  if (scanf("%d%d", &n, &m) != 2)
    return 2;
  Graph<double> g(n);
  while (m--) {
    int u, v;
    double w = 1;
    if (scanf("%d%d", &u, &v) != 2 || ! (0 <= u && u < n && 0 <= v && v < n))
      return 2;
    if (use_w)
      if (scanf("%lf", &w) != 1 || ! (0 <= w))
        return 2;
    g.addEdge(u, v, w);
  }

  vector<Vector<double, 2>> pos(n);
  Circle<double, 2> circle(space);
  circle(g, pos);
  (*algo)(g, pos);
  delete algo;
  for (int u = 0; u < n; u++)
    printf("%.2lf %.2lf\n", pos[u][0], pos[u][1]);

  cmdline_parser_free(&args_info);
}
