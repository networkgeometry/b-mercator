#include "../include/generatingSDbipartiteSD_unix.hpp"

int main(int argc , char *argv[])
{
  generatingSDBipartiteSD_t the_graph;
  if(parse_options(argc, argv, the_graph))
  {
    the_graph.load_hidden_variables();
    the_graph.generate_edgelist();
    the_graph.generate_bipartite_edgelist();
  }
  return EXIT_SUCCESS;
}