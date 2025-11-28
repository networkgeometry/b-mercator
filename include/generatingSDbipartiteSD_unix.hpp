#ifndef GENERATINGSDBIPARTITESD_UNIX_HPP_INCLUDED
#define GENERATINGSDBIPARTITESD_UNIX_HPP_INCLUDED

#include <cstdlib>
#include <string>
#include <unistd.h>
#include "generatingSDbipartiteSD.hpp"


void print_usage()
{
  std::string_view message = R"""(
NAME
  generatingSDbipartiteSD -- a program to generate complex networks in the S^D metric space

SYNOPSIS
  generatingSDbipartiteSD [options]
  )""";
  std::cout << message << std::endl;
}

void print_help()
{
  std::string_view help = R"""(
The following options are available:
  
  // Unipartite network
  -b [BETA]        Specifies the value for parameter beta.
  -g [GAMMA]       Exponent of the power-law distribution for hidden degrees.
  -n [SIZE]        Network size.
  -k [MEAN_DEGREE] Mean degree of nodes.
  -l [KAPPAS]      File consisting of the hidden degrees
  
  // Bipartite network
  -f [FEATURES]          Number of features
  -p [BETA_B]            Specifies the value for parameter beta in the bipartite network
  -q [MEAN_DEGREE_NODES] Mean degree of nodes in the bipartite networks
  -e [GAMMA_N]           Exponent of the power-law distribution for hidden degrees of nodes
  -t [GAMMA_F]           Exponent of the power-law distribution for hidden degrees of features 
  
  // Additional options
  -d [DIMENSION]   Specifies model's dimension (S^D) both for unipartite and bipartite networks
  -c [CORRELATION] Correlation between positions of nodes in unipartite and bipartite networks
  -s [SEED]        Program uses a custom seed for the random number generator. Default: EPOCH.
  -v               Outputs the hidden variables (kappa and nodes'positions) used to the generate the network into a file (uses the edgelist's rootname).
  -h               Print this message on screen and exit.
  -o [FILENAME]    Name of the output file (without extension) (default: net)
  -r [FOLDER]      Name of the output folder 
  )""";
  std::cout << help << std::endl;
}

bool parse_options(int argc , char *argv[], generatingSDBipartiteSD_t &the_graph)
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
      // Unipartite parameters
      case 'b':
        the_graph.BETA = std::stod(optarg);
        break;

      case 'g':
        the_graph.GAMMA = std::stod(optarg);
        break;
      
      case 'n':
        the_graph.NB_VERTICES = std::stoi(optarg);
        break;

      case 'k':
        the_graph.MEAN_DEGREE = std::stod(optarg);
        break;

      case 'l':
        the_graph.HIDDEN_VARIABLES_FILENAME = optarg;
        break;

      // Bipartite parameters 
      case 'f':
        the_graph.NB_FEATURES = std::stoi(optarg);
        break;
      
      case 'p':
        the_graph.BIPARTITE_BETA = std::stod(optarg);
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
      
      case 'd':
        the_graph.DIMENSION = std::stoi(optarg);
        break;

      case 'c':
        the_graph.POS_CORRELATION = std::stod(optarg);
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


#endif // GENERATINGSDBIPARTITESD_UNIX_HPP_INCLUDED