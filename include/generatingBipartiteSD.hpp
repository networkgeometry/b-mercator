#ifndef GENERATINGBIPARTITESD_HPP_INCLUDED
#define GENERATINGBIPARTITESD_HPP_INCLUDED

/**
 * Compile: g++ -O3 -std=c++17 -lboost_system -lboost_math_c99 src/generatingBipartiteSD_unix.cpp -o genSDbipartite
 * 
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <boost/math/quadrature/gauss.hpp>
#include <functional>


double integrate_3d(const std::function<double(double, double, double)>& func, 
                    double ax, double bx, double ay, double by, double az, double bz)
{
    using namespace boost::math::quadrature;
    auto integrand_z = [&](double z) {
        auto integrand_y = [&](double y) {
            auto integrand_x = [&](double x) {
                return func(x, y, z);
            };
            gauss<double, 50> integrator_x;
            return integrator_x.integrate(integrand_x, ax, bx);
        };
        gauss<double, 50> integrator_y;
        return integrator_y.integrate(integrand_y, ay, by);
    };
    gauss<double, 50> integrator_z;
    return integrator_z.integrate(integrand_z, az, bz);
}

class generatingBipartiteSD_t
{
  // Flags controlling options.
  public:
    bool OUTPUT_VERTICES_PROPERTIES = false;
    bool WITH_INPUT_KAPPAS = false;
  // Global parameters.
  public:
    // Random number generator seed.
    int SEED = std::time(NULL);

    // Unipartite parameters
    int NB_VERTICES = -1; // Network size
    
    // Bipartite parameters
    int NB_FEATURES = -1; // Number of features
    double BIPARTITE_BETA = -1; // Parameter beta_b
    double MEAN_DEGREE_NODES = -1; // Mean degree of nodes in bipartite network
    double GAMMA_N = -1; // Exponent of the power-law distribution for hidden degrees of nodes
    double GAMMA_F = -1; // Exponent of the power-law distribution for hidden degrees of features
    double BIPARTITE_MU = -1; 
    double POS_CORRELATION = 1; // Correlation between nodes' positions in unipartite and bipartite networks
    
    // Rootname for the output files;
    std::string OUTPUT_ROOTNAME = "";
    // Input hidden variables filename.
    std::string HIDDEN_VARIABLES_FILENAME = "";
    // Name of the output files
    std::string OUTPUT_FILENAME = "net";
    
    int DIMENSION = 1; // Dimension of the model S^D
  
  private:
    const double PI = 3.141592653589793238462643383279502884197;
    const double NUMERICAL_ZERO = 1e-10;
    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform_01;
    std::normal_distribution<double> normal_01;
    // Mapping the numerical ID of vertices to their name.
    std::vector<std::string> Num2NameNodes;
    std::vector<std::string> Num2NameFeatures;
  // Objects related to the graph ensemble.
  private:
    // Bipartite variables
    std::vector<double> kappa_nodes;
    std::vector<double> kappa_features;
    std::vector<double> theta_nodes;
    std::vector<double> theta_features;
    std::vector<std::vector<double>> d_positions_nodes;
    std::vector<std::vector<double>> d_positions_features;
    
  public:
    generatingBipartiteSD_t() {};
    void generate_edgelist(int width = 15);
    void generate_bipartite_edgelist(int width = 15);
  
  private:
    void save_bipartite_properties(std::vector<int>& rdegree_nodes, std::vector<int>&rdegree_features, 
                                   std::vector<double>& edegree_nodes, std::vector<double>& edegree_features, int width);
    std::vector<double> generate_random_d_vector(int dim, double radius);
    double compute_angle_d_vectors(const std::vector<double> &v1, const std::vector<double> &v2);
    inline double compute_radius(int dim, int N) const;
    std::string get_time();
    void validate_input_parameters();
    void generate_powerlaw_distribution(int n, double gamma, double mean_degree, std::vector<double> &gen_kappa, 
                                        std::vector<std::string> &num2name, std::string prefix="v");
    void normalize_and_rescale_vector(std::vector<double> &v, double radius);
};


void generatingBipartiteSD_t::validate_input_parameters() 
{
  // Either n,gamma,mean_degree is specified or HIDDEN_VARIABLES_FILENAME is not empty
  if (!HIDDEN_VARIABLES_FILENAME.empty()) 
  {
    WITH_INPUT_KAPPAS = true;
  } else {
    WITH_INPUT_KAPPAS = false;
    HIDDEN_VARIABLES_FILENAME = ".";
  }
  auto path = std::filesystem::path(HIDDEN_VARIABLES_FILENAME);
  path = std::filesystem::absolute(path);
  HIDDEN_VARIABLES_FILENAME = path.string();
  if (OUTPUT_ROOTNAME.empty())
  {
    OUTPUT_ROOTNAME = path.parent_path().string() + "/" + OUTPUT_FILENAME;
  } else 
  {
    OUTPUT_ROOTNAME = OUTPUT_ROOTNAME + "/" + OUTPUT_FILENAME;
  }
}

void generatingBipartiteSD_t::generate_powerlaw_distribution(int n, double gamma, double mean_degree, 
                                                               std::vector<double> &gen_kappa, std::vector<std::string> &num2name, 
                                                               std::string prefix)
{ 
  const double kappa_0 = (1 - 1.0 / n) / (1 - std::pow(n, (2.0 - gamma) / (gamma - 1.0))) * (gamma - 2) / (gamma - 1) * mean_degree;
  const double base = 1 - 1.0 / n;
  const double power = 1 / (1 - gamma);
  
  for (int i=0; i<n; ++i) { 
    auto random_kappa = kappa_0 * std::pow(1 - uniform_01(engine) * base, power);
    gen_kappa.push_back(random_kappa);
    num2name.push_back(prefix + std::to_string(i));
  }
}

void generatingBipartiteSD_t::generate_bipartite_edgelist(int width)
{
  std::cout << "HERE\n";
  validate_input_parameters();

  const double mean_degree_features = ((double)NB_VERTICES) / NB_FEATURES * MEAN_DEGREE_NODES;
  
  generate_powerlaw_distribution(NB_VERTICES, GAMMA_N, MEAN_DEGREE_NODES, kappa_nodes, Num2NameNodes, "v"); 
  generate_powerlaw_distribution(NB_FEATURES, GAMMA_F, mean_degree_features, kappa_features, Num2NameFeatures, "f");

  const double radius = compute_radius(DIMENSION, NB_VERTICES);

  std::string edgelist_filename = OUTPUT_ROOTNAME + ".bipartite.edge";
  
  std::vector<double> edegree_nodes;
  std::vector<double> edegree_features;
  std::vector<int> rdegree_nodes;
  std::vector<int> rdegree_features;
  
  // Initializes the containers for the expected and real degrees.
  if(OUTPUT_VERTICES_PROPERTIES)
  {
    edegree_nodes.resize(NB_VERTICES, 0);
    edegree_features.resize(NB_FEATURES, 0);
    rdegree_nodes.resize(NB_VERTICES, 0);
    rdegree_features.resize(NB_FEATURES, 0);
  }

  if (BIPARTITE_BETA > DIMENSION) {
    const auto top = BIPARTITE_BETA * std::tgamma(DIMENSION / 2.0) * std::sin(DIMENSION * PI / BIPARTITE_BETA);
    const auto bottom = MEAN_DEGREE_NODES * 2 * std::pow(PI, 1 + DIMENSION / 2.0);
    BIPARTITE_MU =  top / bottom;

    // Compute mu numerically
    const double kappa_0_nodes = *std::min_element(kappa_nodes.begin(), kappa_nodes.end());
    const double kappa_c_nodes = kappa_0_nodes * pow(NB_VERTICES, 1.0 / (GAMMA_N - 1));
    
    const double kappa_0_features = *std::min_element(kappa_features.begin(), kappa_features.end());
    const double kappa_c_features = kappa_0_features * pow(NB_FEATURES, 1.0 / (GAMMA_F - 1));
    
    const double prefactorA1 = (GAMMA_N - 1) * pow(kappa_0_nodes, GAMMA_N - 1) / (1 - pow(kappa_c_nodes / kappa_0_nodes, 1 - GAMMA_N));
    const double prefactorA2 = (GAMMA_F - 1) * pow(kappa_0_features, GAMMA_F - 1) / (1 - pow(kappa_c_features / kappa_0_features, 1 - GAMMA_F));
    const double prefactorB = tgamma((DIMENSION + 1.0) / 2.0) / ((sqrt(PI) * tgamma(DIMENSION / 2.0)));
    const double prefactor = NB_VERTICES * prefactorA1 * prefactorA2 * prefactorB;
    
    while (true) {
      auto f = [&](double kappa_n, double kappa_f, double theta) {
        const double top = pow(kappa_n, -GAMMA_N) * pow(kappa_f, -GAMMA_F) * pow(sin(theta), DIMENSION - 1);
        const double bottom = 1 + pow(radius * theta / pow(BIPARTITE_MU * kappa_n * kappa_f, 1.0 / DIMENSION), BIPARTITE_BETA);
        return top / bottom;
      };
      double computedMeanDegreeFeatures = prefactor * integrate_3d(f, kappa_0_nodes, kappa_c_nodes, kappa_0_features, kappa_c_features, 0, PI);
      
      if (fabs(computedMeanDegreeFeatures - mean_degree_features) < 0.1)
        break;
      
      if (computedMeanDegreeFeatures < mean_degree_features)
        BIPARTITE_MU *= 1.5;
      else
        BIPARTITE_MU /= 2;
    }

  } else if (BIPARTITE_BETA < DIMENSION) {
    const auto top = (DIMENSION - BIPARTITE_BETA) * std::tgamma(DIMENSION / 2.0);
    const auto bottom = 2 * std::pow(PI, (3 * DIMENSION / 2.0) - BIPARTITE_BETA) * MEAN_DEGREE_NODES * std::pow(NB_VERTICES, 1 - BIPARTITE_BETA / DIMENSION);
    const auto last = std::pow((2 * std::pow(PI, (DIMENSION + 1) / 2.0)) / std::tgamma((DIMENSION + 1) / 2.0), 1 - BIPARTITE_BETA / DIMENSION);
    BIPARTITE_MU = top / bottom * last;
    
    // Compute mu numerically
    const double kappa_0_nodes = *std::min_element(kappa_nodes.begin(), kappa_nodes.end());
    const double kappa_c_nodes = kappa_0_nodes * pow(NB_VERTICES, 1.0 / (GAMMA_N - 1));
    
    const double kappa_0_features = *std::min_element(kappa_features.begin(), kappa_features.end());
    const double kappa_c_features = kappa_0_features * pow(NB_FEATURES, 1.0 / (GAMMA_F - 1));
    
    const double prefactorA1 = (GAMMA_N - 1) * pow(kappa_0_nodes, GAMMA_N - 1) / (1 - pow(kappa_c_nodes / kappa_0_nodes, 1 - GAMMA_N));
    const double prefactorA2 = (GAMMA_F - 1) * pow(kappa_0_features, GAMMA_F - 1) / (1 - pow(kappa_c_features / kappa_0_features, 1 - GAMMA_F));
    const double prefactorB = tgamma((DIMENSION + 1.0) / 2.0) / ((sqrt(PI) * tgamma(DIMENSION / 2.0)));
    const double prefactor = NB_VERTICES * prefactorA1 * prefactorA2 * prefactorB;
    
    while (true) {
      auto f = [&](double kappa_n, double kappa_f, double theta) {
        const double top = pow(kappa_n, -GAMMA_N) * pow(kappa_f, -GAMMA_F) * pow(sin(theta), DIMENSION - 1);
        const double bottom = 1 + pow(radius * theta, BIPARTITE_BETA) / (BIPARTITE_MU * kappa_n * kappa_f);
        return top / bottom;
      };
      double computedMeanDegreeFeatures = prefactor * integrate_3d(f, kappa_0_nodes, kappa_c_nodes, kappa_0_features, kappa_c_features, 0, PI);
      
      std::cout << "computedMeanDegreeFeatures = " << computedMeanDegreeFeatures << " mean_degree_features = " << mean_degree_features << "\n";
      if (fabs(computedMeanDegreeFeatures - mean_degree_features) < 0.1)
        break;
      
      if (computedMeanDegreeFeatures < mean_degree_features)
        BIPARTITE_MU *= 1.5;
      else
        BIPARTITE_MU /= 2;
    }

  } else {
    throw std::invalid_argument("Case beta_b = dimension is not implemented yet.");
  }

  if (POS_CORRELATION == 0)
    POS_CORRELATION = 1e-8; // to avoid nan in tan()

  if (DIMENSION == 1) {
    theta_nodes.clear();
    theta_nodes.resize(NB_VERTICES);

    for (int f=0; f<NB_VERTICES; ++f)
      theta_nodes[f] = 2 * PI * uniform_01(engine);
    
    theta_features.clear();
    theta_features.resize(NB_FEATURES);
    for (int f=0; f<NB_FEATURES; ++f)
      theta_features[f] = 2 * PI * uniform_01(engine);
  } else {
    d_positions_nodes.clear();
    d_positions_nodes.resize(NB_VERTICES);
    for(int f=0; f<NB_VERTICES; ++f)
      d_positions_nodes[f] = generate_random_d_vector(DIMENSION, radius);

    d_positions_features.clear();
    d_positions_features.resize(NB_FEATURES);
    for(int f=0; f<NB_FEATURES; ++f)
      d_positions_features[f] = generate_random_d_vector(DIMENSION, radius);
  }

  // Opens the stream and terminates if the operation did not succeed.
  std::fstream edgelist_file(edgelist_filename.c_str(), std::fstream::out);
  if(!edgelist_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }

  edgelist_file << "# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=" << std::endl;
  edgelist_file << "# Generated on:           " << get_time()                << std::endl;
  edgelist_file << "# Hidden variables file:  " << HIDDEN_VARIABLES_FILENAME << std::endl;
  edgelist_file << "# Seed:                   " << SEED                      << std::endl;
  edgelist_file << "#"                                                       << std::endl;
  edgelist_file << "# Parameters"                                            << std::endl;
  edgelist_file << "#   - dimension:             " << DIMENSION                 << std::endl;
  edgelist_file << "#   - nb. vertices:          " << NB_VERTICES               << std::endl;
  edgelist_file << "#   - nb. features:          " << NB_FEATURES               << std::endl;
  edgelist_file << "#   - beta_b:                " << BIPARTITE_BETA            << std::endl;
  edgelist_file << "#   - mu_b:                  " << BIPARTITE_MU              << std::endl;
  edgelist_file << "#   - radius:                " << radius                    << std::endl;
  edgelist_file << "#   - gamma_n:               " << GAMMA_N                   << std::endl;
  edgelist_file << "#   - gamma_f:               " << GAMMA_F                   << std::endl;
  edgelist_file << "#   - mean degree nodes:     " << MEAN_DEGREE_NODES         << std::endl;
  edgelist_file << "#   - mean degree features:  " << mean_degree_features      << std::endl;
  edgelist_file << "# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=" << std::endl;
  edgelist_file << "#";
  edgelist_file << std::setw(width - 1) << "Vertex" << " ";
  edgelist_file << std::setw(width)     << "Feature" << " ";
  edgelist_file << std::endl;

  double dtheta;
  for(int v=0; v<NB_VERTICES; ++v) {
    for(int f=0; f<NB_FEATURES; ++f) {
      if (DIMENSION == 1) {
        dtheta = PI - std::fabs(PI - std::fabs(theta_nodes[v] - theta_features[f]));
      } else {
        dtheta = compute_angle_d_vectors(d_positions_nodes[v], d_positions_features[f]);
      }
      
      const auto inside_top = std::pow(radius * dtheta, BIPARTITE_BETA);
      const auto inside_bottom = std::pow(std::pow(BIPARTITE_MU * kappa_nodes[v] * kappa_features[f], 1.0 / DIMENSION), std::max((double)DIMENSION, BIPARTITE_BETA));
      const auto prob = 1 / (1 + inside_top / inside_bottom);

      if(uniform_01(engine) < prob)
      {
        edgelist_file << std::setw(width) << Num2NameNodes[v] << " ";
        edgelist_file << std::setw(width) << Num2NameFeatures[f] << " ";
        edgelist_file << std::endl;
        if(OUTPUT_VERTICES_PROPERTIES)
        {
          rdegree_nodes[v] += 1;
          rdegree_features[f] += 1;
        }
      }
      if(OUTPUT_VERTICES_PROPERTIES)
      {
        edegree_nodes[v] += prob;
        edegree_features[f] += prob;
      }
    }
  }
  edgelist_file.close();
  
  if(OUTPUT_VERTICES_PROPERTIES)
  {
    save_bipartite_properties(rdegree_nodes, rdegree_features, edegree_nodes, edegree_features, width);
  }
}


void generatingBipartiteSD_t::save_bipartite_properties(std::vector<int>& rdegree_nodes, std::vector<int>&rdegree_features, 
                                                          std::vector<double>& edegree_nodes, std::vector<double>& edegree_features, int width)
{
  double kappa_min_nodes = *std::min_element(kappa_nodes.begin(), kappa_nodes.end());
  double kappa_min_features = *std::min_element(kappa_features.begin(), kappa_features.end());
  const auto R = compute_radius(DIMENSION, NB_VERTICES);
  const double hyp_radius = 2 * std::log(2 * R / std::pow(BIPARTITE_MU * kappa_min_nodes * kappa_min_features, 1.0 / DIMENSION));

  // Save nodes from the bipartite network
  std::string hidden_variables_nodes_filename = OUTPUT_ROOTNAME + ".gen_coord_bipartite_nodes";
  std::fstream hidden_variables_nodes_file(hidden_variables_nodes_filename.c_str(), std::fstream::out);
  if(!hidden_variables_nodes_file.is_open() )
  {
    std::cerr << "Could not open file: " << hidden_variables_nodes_filename << "." << std::endl;
    std::terminate();
  }

  hidden_variables_nodes_file << "#";
  hidden_variables_nodes_file << std::setw(width - 1) << "Vertex"   << " ";
  hidden_variables_nodes_file << std::setw(width)     << "Kappa"    << " ";
  hidden_variables_nodes_file << std::setw(width)     << "Hyp.Rad."    << " ";
  if (DIMENSION == 1) 
  {
    hidden_variables_nodes_file << std::setw(width)     << "Theta"    << " ";
  } else {
    for (int i=0; i<DIMENSION+1; ++i)
      hidden_variables_nodes_file << std::setw(width)   << "Pos." << i << " ";  
  }
  hidden_variables_nodes_file << std::setw(width)     << "RealDeg." << " ";
  hidden_variables_nodes_file << std::setw(width)     << "Exp.Deg." << " ";
  hidden_variables_nodes_file << std::endl;
  // Writes the hidden variables.
  for(int v(0); v<NB_VERTICES; ++v)
  {
    hidden_variables_nodes_file << std::setw(width) << Num2NameNodes[v]                                                                 << " ";
    hidden_variables_nodes_file << std::setw(width) << kappa_nodes[v]                                                              << " ";
    hidden_variables_nodes_file << std::setw(width) << hyp_radius - (2.0 / DIMENSION) * std::log(kappa_nodes[v] / kappa_min_nodes) << " ";
    if (DIMENSION == 1) {
      hidden_variables_nodes_file << std::setw(width) << theta_nodes[v]                                                            << " ";
    } else {
      for (int i=0; i<DIMENSION+1; ++i)
        hidden_variables_nodes_file << std::setw(width)   << d_positions_nodes[v][i] << " ";
    }
    hidden_variables_nodes_file << std::setw(width) << rdegree_nodes[v]                                                     << " ";
    hidden_variables_nodes_file << std::setw(width) << edegree_nodes[v]                                                     << " ";
    hidden_variables_nodes_file << std::endl;
  }
  hidden_variables_nodes_file.close();


  // Save features from the bipartite network
  std::string hidden_variables_features_filename = OUTPUT_ROOTNAME + ".gen_coord_bipartite_features";
  std::fstream hidden_variables_features_file(hidden_variables_features_filename.c_str(), std::fstream::out);
  if(!hidden_variables_features_file.is_open() )
  {
    std::cerr << "Could not open file: " << hidden_variables_features_filename << "." << std::endl;
    std::terminate();
  }

  hidden_variables_features_file << "#";
  hidden_variables_features_file << std::setw(width - 1) << "Vertex"   << " ";
  hidden_variables_features_file << std::setw(width)     << "Kappa"    << " ";
  hidden_variables_features_file << std::setw(width)     << "Hyp.Rad."    << " ";
  if (DIMENSION == 1) 
  {
    hidden_variables_features_file << std::setw(width)     << "Theta"    << " ";
  } else {
    for (int i=0; i<DIMENSION+1; ++i)
      hidden_variables_features_file << std::setw(width)   << "Pos." << i << " ";  
  }
  hidden_variables_features_file << std::setw(width)     << "RealDeg." << " ";
  hidden_variables_features_file << std::setw(width)     << "Exp.Deg." << " ";
  hidden_variables_features_file << std::endl;
  // Writes the hidden variables.
  for(int f=0; f<NB_FEATURES; ++f)
  {
    hidden_variables_features_file << std::setw(width) << Num2NameFeatures[f]                                                               << " ";
    hidden_variables_features_file << std::setw(width) << kappa_features[f]                                                                 << " ";
    hidden_variables_features_file << std::setw(width) << hyp_radius - (2.0 / DIMENSION) * std::log(kappa_features[f] / kappa_min_features) << " ";
    if (DIMENSION == 1) {
      hidden_variables_features_file << std::setw(width) << theta_features[f]                                                            << " ";
    } else {
      for (int i=0; i<DIMENSION+1; ++i)
        hidden_variables_features_file << std::setw(width)   << d_positions_features[f][i] << " ";
    }
    hidden_variables_features_file << std::setw(width) << rdegree_features[f]                                                     << " ";
    hidden_variables_features_file << std::setw(width) << edegree_features[f]                                                     << " ";
    hidden_variables_features_file << std::endl;
  }
  hidden_variables_features_file.close();

}

std::string generatingBipartiteSD_t::get_time()
{
  // Gets the current date/time.
  time_t theTime = time(NULL);
  struct tm *aTime = gmtime(&theTime);
  int year    = aTime->tm_year + 1900;
  int month   = aTime->tm_mon + 1;
  int day     = aTime->tm_mday;
  int hours   = aTime->tm_hour;
  int minutes = aTime->tm_min;
  // Format the string.
  std::string the_time = std::to_string(year) + "/";
  if(month < 10)
    the_time += "0";
  the_time += std::to_string(month) + "/";
  if(day < 10)
    the_time += "0";
  the_time += std::to_string(day) + " " + std::to_string(hours) + ":";
  if(minutes < 10)
    the_time += "0";
  the_time += std::to_string(minutes) + " UTC";
  // Returns the date/time.
  return the_time;
}

std::vector<double> generatingBipartiteSD_t::generate_random_d_vector(int dim, double radius) {
  std::vector<double> positions;
  positions.resize(dim + 1);
  double norm{0};
  for (auto &pos : positions) {
    pos = normal_01(engine);
    norm += pos * pos;
  }
  norm /= std::sqrt(norm);
  // Normalize vector
  for (auto &pos: positions)
    pos /= norm;

  for (auto &pos: positions)
    pos *= radius;
  return positions;
}

double generatingBipartiteSD_t::compute_angle_d_vectors(const std::vector<double> &v1, const std::vector<double> &v2) {
  double angle{0}, norm1{0}, norm2{0};
  for (int i = 0; i < v1.size(); ++i) {
    angle += v1[i] * v2[i];
    norm1 += v1[i] * v1[i];
    norm2 += v2[i] * v2[i];
  }
  norm1 /= sqrt(norm1);
  norm2 /= sqrt(norm2);
  
  const auto result = angle / (norm1 * norm2);
  if (std::fabs(result - 1) < NUMERICAL_ZERO)
    return 0; // the same vectors
  else
    return std::acos(result);
}

inline double generatingBipartiteSD_t::compute_radius(int dim, int N) const
{
  const auto inside = N / (2 * std::pow(PI, (dim + 1) / 2.0)) * std::tgamma((dim + 1) / 2.0);
  return std::pow(inside, 1.0 / dim);
}

void generatingBipartiteSD_t::normalize_and_rescale_vector(std::vector<double> &v, double radius) {
  int dim = v.size() - 1;
  double norm=0;
  for (int i=0; i<dim + 1; ++i)
    norm += v[i] * v[i];
  
  norm = std::sqrt(norm);
  for (int i=0; i<dim + 1; ++i)
    v[i] /= norm;
  
  for (int i=0; i<dim + 1; ++i)
    v[i] *= radius;
}


#endif // GENERATINGBIPARTITESD_HPP_INCLUDED