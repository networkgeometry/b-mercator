#include "../include/generatingBipartiteSD_unix.hpp"

int main(int argc , char *argv[])
{
  generatingBipartiteSD_t the_graph;
  if(parse_options(argc, argv, the_graph))
  {
    the_graph.generate_bipartite_edgelist();
  }
  return EXIT_SUCCESS;
}