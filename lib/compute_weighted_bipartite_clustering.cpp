#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

int nb_feature_vertices = 0;
int nb_vertices = 0;
std::vector<std::set<int>> bipartite_adjacency_list_nodes;
std::vector<std::set<int>> bipartite_adjacency_list_features;
std::map<std::string, int> Name2Num;
std::map<std::string, int> Name2NumFeature;


void load_bipartite_edgelist(std::string bipartite_edgelist_filename)
{
  std::ifstream bipartite_edgelist_file;
  std::stringstream one_line;
  int v1, v2;
  std::string full_line, name1_str, name2_str;
  std::map< std::string, int >::iterator name_it;
  
  bipartite_adjacency_list_nodes.clear();
  bipartite_adjacency_list_features.clear();

  bipartite_edgelist_file.open(bipartite_edgelist_filename.c_str(), std::ios_base::in);
  if(!bipartite_edgelist_file.is_open())
  {
    std::cerr << "Could not open file: " << bipartite_edgelist_filename<< "." << std::endl;
    std::terminate();
  }
  else
  {
    while (!bipartite_edgelist_file.eof())
    {
      // Reads a line of the file.
      std::getline(bipartite_edgelist_file, full_line); bipartite_edgelist_file >> std::ws;
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      // Skips lines of comment.
      if(name1_str == "#")
      {
        one_line.clear();
        continue;
      }
      one_line >> name2_str >> std::ws;
      one_line.clear();

      if(name1_str != name2_str)
      {
        name_it = Name2Num.find(name1_str);
        if(name_it == Name2Num.end())
        {
          // New vertex
          v1 = nb_vertices;
          Name2Num[name1_str] = v1;
          bipartite_adjacency_list_nodes.emplace_back();
          ++nb_vertices;
        }
        else
        {
          v1 = name_it->second;
        }

        // Is name2 new?
        name_it = Name2NumFeature.find(name2_str);
        if(name_it == Name2NumFeature.end())
        {
          // New feature.
          v2 = nb_feature_vertices;
          Name2NumFeature[name2_str] = v2;
          bipartite_adjacency_list_features.emplace_back();
          ++nb_feature_vertices;
        }
        else
        {
          v2 = name_it->second;
        }
        
        bipartite_adjacency_list_nodes[v1].insert(v2);
        bipartite_adjacency_list_features[v2].insert(v1);
      } 
    }
  }
  bipartite_edgelist_file.close();
}


double compute_weighted_bipartite_clustering(const std::vector<std::set<int>> &bipartite_adjacency_list1, 
                                             const std::vector<std::set<int>> &bipartite_adjacency_list2)
{
  double bipartite_clustering = 0;
  const int size = bipartite_adjacency_list1.size();
  int size_degree_gr2 = size;
#pragma omp parallel for default(shared) reduction(+:bipartite_clustering)
  for (int i=0; i<size; ++i) // type 1 nodes
  {
    double num_triangles = 0;
    const auto neighbours_i = bipartite_adjacency_list1[i]; // neighbours of node i (type 2 nodes) 
    const int neighbours_i_size = neighbours_i.size();
    if (neighbours_i_size < 2) {
      size_degree_gr2 -= 1;
      continue; 
    }
    
    for (const int n1: neighbours_i)
    {
      for (const int n2: neighbours_i) 
      {
        if (n1 == n2)
          continue;
        
        const auto neighbours_j = bipartite_adjacency_list2[n1]; // neighbours of neighbour of node i (type 1 nodes)
        const auto neighbours_k = bipartite_adjacency_list2[n2];
        // Compute the common type 1 nodes between two sets
        std::set<int> intersect;
        std::set_intersection(neighbours_j.begin(), neighbours_j.end(), 
                              neighbours_k.begin(), neighbours_k.end(),
                              std::inserter(intersect, intersect.begin()));
        int common_nodes = intersect.size();
        if (common_nodes > 1) // ignore the node i from the computations
          num_triangles += (common_nodes - 1) / (double)(std::min(neighbours_j.size(), neighbours_k.size()) - 1) / 2; // do not count twice
      }
    }
    const double tmp_clustering = (2 * num_triangles) / ((neighbours_i_size) * (neighbours_i_size - 1));
    bipartite_clustering += tmp_clustering;
  }
  return bipartite_clustering / size_degree_gr2; 
}


int main(int argc , char *argv[]) {
  const std::string bipartite_edgelist_path = argv[1];

  load_bipartite_edgelist(bipartite_edgelist_path);
  
  const double cb_nodes = compute_weighted_bipartite_clustering(bipartite_adjacency_list_nodes, bipartite_adjacency_list_features);
  const double cb_features = compute_weighted_bipartite_clustering(bipartite_adjacency_list_features, bipartite_adjacency_list_nodes);

  std::cout << cb_nodes << "," << cb_features << std::endl;
}