#ifndef GENERATINGBIPARTITESD_UNIX_HPP_INCLUDED
#define GENERATINGBIPARTITESD_UNIX_HPP_INCLUDED

#include <cstdlib>
#include <string>
#include <unistd.h>
#include "generatingBipartiteSD.hpp"


void print_usage()
{
  std::string_view message = R"""(
NAME
  generatingBipartiteSD -- a program to generate synthetic networks from the bipartite-S^D model

SYNOPSIS
  generatingBipartiteSD [options]
  )""";
  std::cout << message << std::endl;
}

void print_help()
{
  std::string_view help = R"""(
The following options are available:
  
  -d [DIMENSION]         Specifies model's dimension (S^D) both for unipartite and bipartite networks
  -n [SIZE]              Number of nodes.
  -f [FEATURES]          Number of features
  -b [BETA_B]            Specifies the value for parameter beta in the bipartite network
  -q [MEAN_DEGREE_NODES] Mean degree of nodes in the bipartite networks
  -e [GAMMA_N]           Exponent of the power-law distribution for hidden degrees of nodes
  -t [GAMMA_F]           Exponent of the power-law distribution for hidden degrees of features 
  
  // Additional options
  -s [SEED]        Program uses a custom seed for the random number generator. Default: EPOCH.
  -v               Outputs the hidden variables (kappa and nodes'positions) used to the generate the network into a file (uses the edgelist's rootname).
  -h               Print this message on screen and exit.
  -o [FILENAME]    Name of the output file (without extension) (default: net)
  -r [FOLDER]      Name of the output folder 
  )""";
  std::cout << help << std::endl;
}

bool parse_options(int argc , char *argv[], generatingBipartiteSD_t &the_graph)
{
  if(argc == 1)
  {
    print_usage();
    print_help();
    return false;
  }

  // Parsing options.
  int opt;
  while ((opt = getopt(argc,argv,"b:d:g:n:k:l:f:p:q:e:t:hs:vo:r:c:")) != -1)
  {
    switch(opt)
    {

      case 'd':
        the_graph.DIMENSION = std::stoi(optarg);
        break;

      case 'b':
        the_graph.BIPARTITE_BETA = std::stod(optarg);
        break;
      
      case 'n':
        the_graph.NB_VERTICES = std::stoi(optarg);
        break;

      case 'f':
        the_graph.NB_FEATURES = std::stoi(optarg);
        break;
      
      case 'q':
        the_graph.MEAN_DEGREE_NODES = std::stod(optarg);
        break;
      
      case 'e':
        the_graph.GAMMA_N = std::stod(optarg);
        break;
     
      case 't':
        the_graph.GAMMA_F = std::stod(optarg);
        break;
      
      // Additional parameters
      case 's':
        the_graph.SEED = std::stoi(optarg);
        break;

      case 'v':
        the_graph.OUTPUT_VERTICES_PROPERTIES = true;
        break;
      
      case 'o':
        the_graph.OUTPUT_FILENAME = optarg;
        break;

      case 'r':
        the_graph.OUTPUT_ROOTNAME = optarg;
        break;

      case 'h':
        print_usage();
        print_help();
        return false;

      default:
        print_usage();
        print_help();
        return false;
    }
  }

  if (the_graph.OUTPUT_ROOTNAME.empty()) {
    std::cout << "the_graph.OUTPUT_ROOTNAME = " << the_graph.OUTPUT_ROOTNAME << std::endl;
    // Uses the default rootname for output files.
    size_t lastdot = the_graph.HIDDEN_VARIABLES_FILENAME.find_last_of(".");
    if(lastdot == std::string::npos)
    {
      the_graph.OUTPUT_ROOTNAME = the_graph.HIDDEN_VARIABLES_FILENAME;
    }
    the_graph.OUTPUT_ROOTNAME = the_graph.HIDDEN_VARIABLES_FILENAME.substr(0, lastdot);
  }

  return true;
}


#endif // GENERATINGBIPARTITESD_UNIX_HPP_INCLUDED