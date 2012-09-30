#include "Core.hh"
#include "FruchtermanReingold.hh"
#include "Cmdline.h"

int main(int argc, char* argv[])
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0)
    return 1;

  ForceDirectedDrawing<double, 2>* algo = NULL;
  switch (args_info.algorithm_arg) {
  case 0:
    {
      auto a = new FruchtermanReingold<double, 2>();
      a->space[0] = args_info.width_arg;
      a->space[1] = args_info.height_arg;
      algo = a;
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
    if (scanf("%d%d", &u, &v) != 2 || ! (0 <= u && u < n && 0 <= v && v < n))
      return 2;
    g.addEdge(u, v, 1);
  }

  vector<Vector<double, 2>> pos = (*algo)(g);
  delete algo;
  for (int u = 0; u < n; u++)
    printf("%.2lf %.2lf\n", pos[u][0], pos[u][1]);

  cmdline_parser_free(&args_info);
}
