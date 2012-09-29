#include "Core.hh"
#include "FruchtermanReingold.hh"
#include "Cmdline.h"

int main(int argc, char* argv[])
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0)
    return 1;

  ForceDirectedAlgorithm* algo = NULL;
  switch (args_info.algorithm_arg) {
  case 0:
    algo = new FruchtermanReingold();
    break;
  default:
    return 2;
  }
  algo->setSize(args_info.width_arg, args_info.height_arg);

  int n, m;
  if (scanf("%d%d", &n, &m) != 2)
    return 2;
  Graph g(n);
  while (m--) {
    int u, v;
    if (scanf("%d%d", &u, &v) != 2 || ! (0 <= u && u < n && 0 <= v && v < n))
      return 2;
    g.addEdge(u, v);
  }

  vector<Vector> pos = (*algo)(g);
  delete algo;
  for (int u = 0; u < n; u++)
    printf("%.2lf %.2lf\n", pos[u].x, pos[u].y);

  cmdline_parser_free(&args_info);
}
