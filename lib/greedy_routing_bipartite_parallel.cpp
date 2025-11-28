#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <complex>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <filesystem>
#include <list>
#include <queue>


using AdjacencyList = std::map<std::string, std::vector<std::string>>;

AdjacencyList bipartite_adjacency_list_nodes;
AdjacencyList bipartite_adjacency_list_features;

const std::string suffixNodes = "NODES";
const std::string suffixFeatures = "FEATURES";


template<typename K, typename V>
std::vector<K> getKeys(const std::map<K, V>& map) {
    std::vector<K> keys;
    keys.reserve(map.size());
    for (const auto& pair : map) {
        keys.push_back(pair.first);
    }
    return keys;
}


template<typename K>
K getRandomKey(const std::vector<K>& keys) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, keys.size() - 1);
    int randomIndex = dis(gen);
    return keys[randomIndex];
}

void print_help() {
    constexpr auto help_message = R"(
        ./greedy_routing_bipartite [dim] [path_to_coords_nodes] [path_to_coords_features] [path_to_bipartite_edgelist] [n_runs]
        
        dim                                   -- dimension of S^D model
        path_to_coords_nodes                  -- path to the file with the nodes' coordinates (.inf_coord_nodes)
        path_to_coords_features               -- path to the file with the features' coordinates (.inf_coord_features)
        path_to_bipartite_edgelist            -- path to the bipartite edgelist (.edge) 
        n_runs                                -- number of greedy routing rounds (default=10*network size)
    )";
    std::cout << help_message << std::endl;
}

void load_coords(int dim,
                 const std::string &coords_path,
                 std::map<std::string, double> &radii,
                 std::map<std::string, double> &thetas, 
                 std::map<std::string, std::vector<double>> &positions,
                 const std::string suffix) {

    std::stringstream one_line;
    std::string full_line, name1_str, name2_str, name3_str;

    std::fstream hidden_variables_file(coords_path.c_str(), std::fstream::in);
    if( !hidden_variables_file.is_open() )
    {
        std::cerr << "Could not open file: " << coords_path << "." << std::endl;
        std::terminate();
    }

    while(!hidden_variables_file.eof()) {
        // Reads a line of the file.
        std::getline(hidden_variables_file, full_line);
        hidden_variables_file >> std::ws;
        one_line.str(full_line);
        one_line >> std::ws;
        one_line >> name1_str >> std::ws;
        // Skips lines of comment.
        if(name1_str == "#")
        {
            one_line.clear();
            continue;
        }

        one_line >> name2_str >> std::ws; // omit kappas

        std::string entity_name = suffix + name1_str;

        if (dim == 1) {
            one_line >> name3_str >> std::ws;
            thetas[entity_name] = std::stod(name3_str);
            one_line >> name3_str >> std::ws;
            if (name3_str == "inf") {// omit zero degree nodes
                thetas.erase(entity_name);
                one_line.clear();
                continue;
            }
            radii[entity_name] = std::stod(name3_str);
        } else {
            one_line >> name3_str >> std::ws;
            if (name3_str == "inf") {// omit zero degree nodes
                one_line.clear();
                continue;
            }
            radii[entity_name] = std::stod(name3_str);

            std::vector<double> tmp_position;
            for (int i=0; i<dim+1; ++i) {
                one_line >> name3_str >> std::ws;
                tmp_position.push_back(std::stod(name3_str));
            }
            positions[entity_name] = tmp_position;
        }
        one_line.clear();
    }
    hidden_variables_file.close();
}


void load_bipartite_edgelist(const std::string& bipartite_edgelist_path, 
                             const std::vector<std::string> &nodesList,
                             const std::vector<std::string> &featuresList) {
    std::stringstream one_line;
    std::string full_line, source_str, target_str;

    std::fstream edgelist_file(bipartite_edgelist_path.c_str(), std::fstream::in);
    if( !edgelist_file.is_open() )
    {
        std::cerr << "Could not open file: " << bipartite_edgelist_path << "." << std::endl;
        std::terminate();
    }

    while(!edgelist_file.eof()) {
        std::getline(edgelist_file, full_line);
        edgelist_file >> std::ws;
        std::istringstream ss(full_line);
        ss >> source_str >> target_str;

        if(source_str == "#" || target_str == "#")
            continue; 

        std::string nodeName = suffixNodes + source_str;
        std::string featureName = suffixFeatures + target_str;

        bool node_found = (std::find(nodesList.begin(), nodesList.end(), nodeName) != nodesList.end());
        bool feature_found = (std::find(featuresList.begin(), featuresList.end(), featureName) != featuresList.end());

        if (node_found && feature_found) {
            bipartite_adjacency_list_nodes[nodeName].push_back(featureName);
            bipartite_adjacency_list_features[featureName].push_back(nodeName);
        }
    }
}

double compute_angle_d_vectors(const std::vector<double> &v1, const std::vector<double> &v2) {
  double angle{0}, norm1{0}, norm2{0};
  for (int i = 0; i < v1.size(); ++i) {
    angle += v1[i] * v2[i];
    norm1 += v1[i] * v1[i];
    norm2 += v2[i] * v2[i];
  }
  norm1 /= sqrt(norm1);
  norm2 /= sqrt(norm2);
  
  const auto result = angle / (norm1 * norm2);
  if (std::fabs(result - 1) < 1e-15)
    return 0;
  else
    return std::acos(result);
}

double compute_distance(int dim, 
                        const std::string &v1, 
                        const std::string &v2, 
                        const std::map<std::string, double> &radiiA, 
                        const std::map<std::string, double> &radiiB,
                        const std::map<std::string, double> &thetasA, 
                        const std::map<std::string, double> &thetasB,
                        const std::map<std::string, std::vector<double>> &positionsA,
                        const std::map<std::string, std::vector<double>> &positionsB) {
    
    double delta_theta = 0;
    if (dim == 1)
        delta_theta = M_PI - std::fabs(M_PI - std::fabs(thetasA.at(v1) - thetasB.at(v2)));
    else
        delta_theta = compute_angle_d_vectors(positionsA.at(v1), positionsB.at(v2));

    double r1 = radiiA.at(v1), r2 = radiiB.at(v2);
    if (delta_theta == 0) {
        return std::fabs(r1 - r2); 
    } else {
        auto dist = 0.5 * ((1 - std::cos(delta_theta)) * std::cosh(r1 + r2) + (1 + std::cos(delta_theta)) * std::cosh(r1 - r2));
        return std::acosh(dist);
    }
}


int bfs_shortest_path(const std::string& start, const std::string& goal) {
    std::queue<std::pair<std::string, int>> to_visit;
    std::set<std::string> visited;

    to_visit.push({start, 0});
    visited.insert(start);

    while (!to_visit.empty()) {
        auto current = to_visit.front();
        to_visit.pop();

        const std::string& node = current.first;
        int distance = current.second;

        if (node == goal) {
            return distance;
        }

        const auto& neighbors = node.find(suffixNodes) != std::string::npos ? 
                                bipartite_adjacency_list_nodes[node] : bipartite_adjacency_list_features[node];

        for (const auto& neighbor : neighbors) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                to_visit.push({neighbor, distance + 1});

                // Add the neighbors of the neighbor, enabling node-to-node and feature-to-feature transitions
                const auto& next_neighbors = neighbor.find(suffixNodes) != std::string::npos ? 
                                             bipartite_adjacency_list_features[neighbor] : bipartite_adjacency_list_nodes[neighbor];
                for (const auto& next_neighbor : next_neighbors) {
                    if (visited.find(next_neighbor) == visited.end()) {
                        visited.insert(next_neighbor);
                        to_visit.push({next_neighbor, distance + 2});
                    }
                }
            }
        }
    }
    return -1; // No path found
}


std::string run_modified_bipartite_greedy_routing(int dim,
                                                  const std::vector<std::string>& nodesA, 
                                                  const std::vector<std::string>& nodesB,
                                                  const std::map<std::string, double>& radii_nodes,
                                                  const std::map<std::string, double>& thetas_nodes,
                                                  const std::map<std::string, std::vector<double>>& positions_nodes,
                                                  const std::map<std::string, double>& radii_features,
                                                  const std::map<std::string, double>& thetas_features,
                                                  const std::map<std::string, std::vector<double>>& positions_features,
                                                  int n_runs) {
    double p_s = 0;
    double mean_strech = 0;
    double max_strech = 0;
    double mean_hop_length = 0;

#pragma omp parallel for reduction(+:p_s,mean_strech,mean_hop_length)
    for (int i=0; i<n_runs; ++i) {
        std::string source = "";
        std::string target = "";
        while (source == target) {
            source = getRandomKey(nodesA);
            target = getRandomKey(nodesB);
        }

        const std::string original_source = source;
        bool is_target_node = target.find(suffixNodes) != std::string::npos;

        std::vector<std::string> hops = {source};
        bool is_package_dropped = false;
        while (source != target) {
            double smallest_distance = 999999999;
            std::string new_source = source;

            // Check kind of nodes A or B 
            bool is_source_node = source.find(suffixNodes) != std::string::npos;
            const auto& neighbors = is_source_node ? bipartite_adjacency_list_nodes.at(source) : bipartite_adjacency_list_features.at(source);

            for (const auto &n: neighbors) {

                if (n == target) {
                    hops.push_back(target);
                    goto found_target2;
                }

                bool is_n_node = n.find(suffixNodes) != std::string::npos; 
                double distance_n_t = 0;
                
                if (!((!is_target_node) != (!is_n_node))) {
                    if (is_target_node) // target - node, n - node
                        distance_n_t = compute_distance(dim, n, target, radii_nodes, radii_nodes, thetas_nodes, thetas_nodes, positions_nodes, positions_nodes);
                    else // target - feature, n - feature
                        distance_n_t = compute_distance(dim, n, target, radii_features, radii_features, thetas_features, thetas_features, positions_features, positions_features);
                } else {
                    if (is_target_node) // target - node, n - feature
                        distance_n_t = compute_distance(dim, n, target, radii_features, radii_nodes, thetas_features, thetas_nodes, positions_features, positions_nodes);
                    else // target - feature, n - node
                        distance_n_t = compute_distance(dim, target, n, radii_features, radii_nodes, thetas_features, thetas_nodes, positions_features, positions_nodes);
                }

                if (distance_n_t < smallest_distance) {
                    new_source = n;
                    smallest_distance = distance_n_t;
                }
            }
            source = new_source;

            if (hops.size() > 1) {
                if (hops[hops.size() - 2] == source) {
                    is_package_dropped = true;
                    break;
                }
            }
            hops.push_back(source);
        }
        found_target2:

        if (!is_package_dropped) {
            ++p_s;
            mean_hop_length += (double)hops.size();
            double stretch = (double)hops.size() / bfs_shortest_path(original_source, target);
            mean_strech += stretch;
            #pragma omp critical
            {
                if (stretch > max_strech)
                    max_strech = stretch;
            }
        }
    }
    mean_strech /= p_s;
    mean_hop_length /= p_s;
    p_s /= n_runs;
    return std::to_string(p_s) + "," + std::to_string(mean_hop_length) + "," + std::to_string(mean_strech) + "," + std::to_string(max_strech);
}


int main(int argc , char *argv[]) {
    if (argc < 5) {
        std::cout << "Error. Wrong number of parameters." << std::endl;
        print_help();
    }

    int dim = std::stoi(argv[1]);
    std::string nodes_coords_path = argv[2];
    std::string features_coords_path = argv[3];
    std::string bipartite_edgelist_path = argv[4];

    std::map<std::string, double> radii_nodes;
    std::map<std::string, double> thetas_nodes;
    std::map<std::string, std::vector<double>> positions_nodes;

    std::map<std::string, double> radii_features;
    std::map<std::string, double> thetas_features;
    std::map<std::string, std::vector<double>> positions_features;

    load_coords(dim, nodes_coords_path, radii_nodes, thetas_nodes, positions_nodes, suffixNodes);
    load_coords(dim, features_coords_path, radii_features, thetas_features, positions_features, suffixFeatures);

    std::vector<std::string> nodesList = getKeys(radii_nodes);
    std::vector<std::string> featuresList = getKeys(radii_features);

    load_bipartite_edgelist(bipartite_edgelist_path, nodesList, featuresList);
    
    // 4 types of GR
    // 1) node    -> node
    // 2) node    -> feature
    // 3) feature -> node
    // 4) feature -> feature

    int n_runs = 10000;
    if (argc == 6)
        n_runs = std::stoi(argv[5]);

    
    // 1) node -> feature
    std::string gr_node_feature = run_modified_bipartite_greedy_routing(dim, nodesList, featuresList, 
                                                                        radii_nodes, thetas_nodes, positions_nodes,
                                                                        radii_features, thetas_features, positions_features, n_runs);

    
    // 2) feature -> node
    std::string gr_feature_node = run_modified_bipartite_greedy_routing(dim, nodesList, featuresList, 
                                                                        radii_nodes, thetas_nodes, positions_nodes,
                                                                        radii_features, thetas_features, positions_features, n_runs);
    
    // 3) node -> node
    std::string gr_node_node = run_modified_bipartite_greedy_routing(dim, nodesList, nodesList, 
                                                                     radii_nodes, thetas_nodes, positions_nodes,
                                                                     radii_features, thetas_features, positions_features, n_runs);

    // 4) feature -> feature
    std::string gr_feature_feature = run_modified_bipartite_greedy_routing(dim, featuresList, featuresList, 
                                                                           radii_nodes, thetas_nodes, positions_nodes,
                                                                           radii_features, thetas_features, positions_features, n_runs);

    std::cout << "type,p_s,mean_hop_length,mean_strech,max_strech" << std::endl;
    std::cout << "node-feature," + gr_node_feature << std::endl;
    std::cout << "feature-node," + gr_feature_node << std::endl;
    std::cout << "node-node," + gr_node_node << std::endl;
    std::cout << "feature-feature," + gr_feature_feature << std::endl; 
}
    