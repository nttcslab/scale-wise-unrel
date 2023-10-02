# Exact and Efficient Network Reliability Evaluation per Outage Scale

This repository includes the codes to reproduce the results of experiments in the paper "Exact and Efficient Network Reliability Evaluation per Outage Scale."

## Requirements

We use [TdZdd](https://github.com/kunisura/TdZdd) for implementing both proposed and baseline methods. We also use [SAPPOROBDD](https://github.com/Shin-ichi-Minato/SAPPOROBDD) for implementing the baseline method. Before building our code, you must place the header files of [TdZdd](https://github.com/kunisura/TdZdd) into this directory. More specifically, all header files of [TdZdd/include/tdzdd](https://github.com/kunisura/TdZdd/tree/master/include/tdzdd) must be placed on `tdzdd/` directory; e.g. `tdzdd/DdEval.hpp`, `tdzdd/DdSpec.hpp`, etc. If you want to build the baseline method, you should build [SAPPOROBDD](https://github.com/Shin-ichi-Minato/SAPPOROBDD); after that, all header files must be placed on `SAPPOROBDD/include/` and the library file `BDD64.a` must be placed on `SAPPOROBDD/lib/`.

After that, if your environment has CMake version >=3.8, you can build all codes with the following commands:

```shell
(moving to src/ directory)
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running this, all binaries are generated in `release/` directory. Instead of `make`, you can build individual binary by the following commands:
```shell
make main  // Building the proposed method; SAPPOROBDD is not needed
make base  // Building the baseline method
```

### Verified environments

We verified that the building process of our codes and the commands presented below worked fine in the following macOS and Linux environments:

- macOS Big Sur 11.2.1 + Apple clang 12.0.0
- CentOS 7.9 + gcc 4.8.5

## How to reproduce experimental results

After building our code, the following binaries are generated: `main` and `base`. `main` implements our proposed method, while `base` implements the baseline method.

All data used in our experiments are in `data.tar.gz`. After extracting, `*.txt` describes the graph files (as an edgelist), and `*.txt.prob` specifies each link's availability. In addition, `src/*.txt.bc<n>.src` specifies the server nodes when $|T|=$`<n>`.

### Experiments

The proposed method can be executed by the following command:

```shell
./main [graph_file] [probability_file] [server_file] [order_file] <client_file>
```

`[graph_file]`, `[probability_file]`, and `[server_file]` specify the path to the file describing the edgelist of the graph, each link's availability, and the list of server nodes, respectively. `[order_file]` specifies the path to the file describing the ordering among links. Note that we already computed the link order by the path-decomposition heuristics and wrote it on `*.txt`, so `[order_file]` can be specified with the same path as `[graph_file]` when reproducing the experimental results. `<client_file>` is an optional argument that specify the path to the file describing the list of client nodes. If `<client_file>` is not specified, clients are set to all the nodes.

After execution, the program computes the probability that exactly $n^\prime$ nodes are disconnected from any server node for every $n^\prime$. Note that this is equivalent to compute $U(n^\prime)$ values in our paper when the client nodes are all the nodes other than servers.

_Running example:_ To compute the scale-wise unreliability of UsCarrier topology when $|T|=5$, run:

```shell
./main ../data/0189-real-UsCarrier.edgelist.txt ../data/0189-real-UsCarrier.edgelist.txt.prob ../data/src/0189-real-UsCarrier.edgelist.txt.bc5.src ../data/src/0189-real-UsCarrier.edgelist.txt
```

The baseline method can also be executed in the same way as follows:

```shell
./base [graph_file] [probability_file] [server_file] [order_file]
```

## License

This software is released under the NTT license, see `LICENSE.pdf`.
