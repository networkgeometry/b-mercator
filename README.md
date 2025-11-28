![B-Mercator](b-mercator-logo.png)

_Mapping bipartite networks into the multidimensional hyperbolic spaces_

## Installation

### Command line executable

Requirements
* A C++17 (or newer) compliant compiler
* `cmake` >= 3.20
* The header `unistd.h`.

```
# Unix (Linux / MAC OS)
chmod +x build.sh
./build.sh -b Release
```

### Docker

See `run_dmercator_docker.py` for more information.

---

## Usage

Given a **bipartite** edgelist `net.edge`, the embedding can be performed as follows
```
./bmercator -d 1 -x -m net.edge
```

where 
* `-d 1` specifies the dimension of the embeddings (here $D=1$, i.e., similarity space is a circle)
* `-x` Bipartite mode. Ignore the unipartite network and only embed the bipartite network.
* `-m [FILENAME]`  Path to the bipartite edgelist


Given a **unipartite** edgelist `net.edge`, the embedding can be performed as follows
```
./bmercator -d 2 net.edge
```
Note that for unipartite network the option `-x` and `-m` are not needed.


### Format of input files

The structure of the network to be embedded is passed to the program via a file containing its edgelist (one link per line). The edgelist file consists in a simple text file with the following convention
```
# lines beginning with "#" are ignored (comments).
# note that nodes' name must be separated by at least one white space.
# there may be white space at the beginning of a line.
[name of node1]  [name of node2]  [remaining information will be ignored]
[name of node2]  [name of node3]  [remaining information will be ignored]
...
```
Note that the nodes' name will be imported as `std::string` and can therefore be virtually anything as long as they do not include white spaces (i.e., there is not need for the nodes to be identified by contiguous integers).

**IMPORTANT**: this class only considers **simple undirected** networks **without self-loops**. Any multiple edges (e.g., if the graph is originally directed) or self-loops will be ignored.

**IMPORTANT**: in the actual version of the code, the network must have **only one component**.


### Output files

For *bipartite* networks, the following files are generated:

* `*.inf_coord_nodes`: Contains information about the embedding procedure, the inferred parameters and the inferred positions of type-A nodes.
* `*.inf_coord_features`: Contains information about the embedding procedure, the inferred parameters and the inferred positions of type-B nodes.
* `*.inf_log`: Contains detailled information about the embedding procedure.

For *unipartite* networks, the following files are generated:

* `*.inf_coord`: Contains information about the embedding procedure, the inferred parameters and the inferred positions.
* `*.inf_log`: Contains detailled information about the embedding procedure.

### Options
```
    -d [DIMENSION] Dimension of the embeddings.
    -m [FILENAME]  Path to the bipartite edgelist for the nodes' features
    -x             Bipartite mode. Ignore the unipartite network and only embed the bipartite
                   network.
    -b [VALUE]     Specify the value for beta to be used for the embedding. By 
                   default the program infers the value of beta based on the average
                   local clustering coefficient of the original edgelist.
    -c             Clean mode. Writes the inferred coordinates in clean file without
                   any lines starting by # to facilitate their use in custom computer
                   programs.
    -f             Fast mode. Does not infer the positions based on likelihood
                   maximization, rather uses only the EigenMap method.
    -k             No post-processing of the values of kappa based on the inferred
                   angular positions (theta) resulting in every vertices with the same
                   degree ending at the same radial position in the hyperbolic disk.
    -o [ROOTNAME]  Specify the rootname used for all output files. Default: uses the
                   rootname of the edgelist file as (i.e., rootname.edge).
    -r [FILENAME]  Refine mode. Reads the inferred positions from a previous run of
                   this program (file *.inf_coord) and refines the inferred positions.
    -q             Quiet mode. Program does not output information about the network
                   and the embedding procedure.
    -s [SEED]      Program uses a custom seed for the random number generator.
                   Default: EPOCH.
    -v             Validation mode. Validates and characterizes the inferred random
                   network ensemble.
```

#### Bipartite mode vs unipartite mode

When embedding a bipartite network, the option `-x` needs to be provided to indicate that the unipartite network should be ignored and only the bipartite network should be embedded. In addition, the path to the edgelist of the bipartite network must be provided using the option `-m`.

```
./bmercator -d <dimension_value> -x -m <bipartite_edgelist_filename>
```

When embedding a unipartite network, neither the option `-x` nor the option `-m` are needed.

```
./bmercator -d <dimension_value> <unipartite_edgelist_filename>
```


#### Dimension
In order to set the dimension of the embedding a parameter `d` need to be set

```
# Command line
./bmercator -d <dimension_value> <edgelist_filename>
```

#### Custom value for beta

A custom value for the parameter `beta` can be provided. By default a value for beta is inferred to reproduce the average local clustering coefficient of the original edgelist.

```
./bmercator -b <beta_value> <edgelist_filename>
```

#### Custom value for the seed of the random number generator

A custom seed for the random number generator can be provided (useful when several embeddings are performed in parallel). By default, `EPOCH` is used.

```
./bmercator -s <seed_value> <edgelist_filename>
```


#### Custom output filename

All generated files are named `<output_rootname>.<extension>` and a custom `<output_rootname>` can be provided. If none is provided, the `<output_rootname>` is extracted from the `<edgelist_filename>` by removing its extension, otherwise the full `<edgelist_filename>` is used as `<output_rootname>` if `<edgelist_filename>`does not have any extension.

```
./bmercator -o <custom_output_rootname> <edgelist_filename>
```


#### Clean output mode

Outputs a file with extension `*.inf_coord_raw` containing the columns 2, 3 and 4 of the file with extension `*.inf_coord`. Rows follow the same order as in the file with extension `*.inf_coord`. The global parameters (i.e., beta, mu, etc.) ate not included in the file. Default is **`false`**.

```
./bmercator -c <edgelist_filename>
```


#### Fast mode

Skip the likelihood maximization step (i.e., only infers the positions using the EigenMap methods). Default is **`false`**. Only applicable where dimension is set to 1.

```
./bmercator -f <edgelist_filename>
```

#### Post-processing of the inferred values of the radial positions

The inferred radial positions are updated based on the inferred angular positions. When deactivated, nodes with the same degree have the same radial position in the hyperbolic disk. Default is **`true`**.

```
./bmercator -k <edgelist_filename>
```


#### Refine mode

When a file containing the previously inferred coordinates is provided (`*.inf_coord`), the program uses the already inferred positions and parameters as a starting point and perform another round of the likelihood maximization step to refine the inferred positions. The use of a different output filename is recommended to keep track of the changes. Default is **`false`**.

```
./bmercator -r <inferred_coordinates_filename> <edgelist_filename>
```

#### Validation mode

Validates and characterizes the inferred random network ensemble. This is done by generating a large number of networks based on the inferred parameters and positions. The following files are generated

* `*.inf_inf_pconn`: Compares the inferred probability of connection with the theoretical one based on the inferred parameters.
* `*.inf_theta_density`: Density the angular distribution
* `*.inf_vprop`: Contains the following topological properties for every nodes in the inferred ensemble and in the original network:
    * degree
    * sum of the degree of the neighbors
    * average degree of the neighbors
    * number of triangles
    * local clustering coefficient
* `*.inf_vstat`: Contains the following statistics of the inferred ensemble.
    * degree distribution
    * spectrum of the average degree of neighbors
    * clustering spectrum
* `*.obs_vstat`: Contains the same statistics as above but for the original network.

```
./bmercator -v <edgelist_filename>
```

See `notebooks/paper-SI-bipartite-validation-synthetic-networks.ipynb` for some plotting examples.


## Publications

- _Mapping bipartite networks into multidimensional hyperbolic spaces_<br>
  [Robert Jankowski](https://robertjankowski.github.io/),
  [Roya Aliakbarisani](https://scholar.google.com/citations?user=99FmtwEAAAAJ&hl=en),
  [M. Ángeles Serrano](http://morfeo.ffn.ub.es/~mariangeles/ws_en/) 
  and [Marián Boguñá](http://complex.ffn.ub.es/~mbogunya/) <br>
  Accepted in Communications Physics <br>
  [Full text]() | [arXiv](https://arxiv.org/abs/2503.04316)

- _The D-Mercator method for the multidimensional hyperbolic embedding of real networks_<br>
  [Robert Jankowski](https://robertjankowski.github.io/),
  [Antoine Allard](http://antoineallard.info),
  [Marián Boguñá](http://complex.ffn.ub.es/~mbogunya/) and 
  [M. Ángeles Serrano](http://morfeo.ffn.ub.es/~mariangeles/ws_en/) <br>
  Nature Commmunications 14, 7585 (2023) <br>
  [Full text](https://www.nature.com/articles/s41467-023-43337-5) | [arXiv](https://arxiv.org/abs/2304.06580)

- _Mercator: uncovering faithful hyperbolic embeddings of complex networks_<br>
[Guillermo García-Pérez*](https://scholar.google.es/citations?user=MibFSJIAAAAJ&hl=en),
[Antoine Allard*](http://antoineallard.info),
[M. Ángeles Serrano](http://morfeo.ffn.ub.es/~mariangeles/ws_en/) and
[Marián Boguñá](http://complex.ffn.ub.es/~mbogunya/)<br>
New Journal of Physics 21, 123033 (2019)<br>
[Full text](https://doi.org/10.1088/1367-2630/ab57d2) | [arXiv](https://arxiv.org/abs/1904.10814)
