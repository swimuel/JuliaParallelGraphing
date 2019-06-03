# SE751-Assignment

before initializing Julia run the following command in the shell to create n threads at startup
1. `export JULIA_NUM_THREADS=n`

installation into Julia from github
1. Open Julia command line
1. Enter the pkg browser by pressing [
1. `add https://github.com/swimuel/SE751-Assignment.git` - as the repo is private you will need to log in to your github account.
1. Exit out of the pkg browser
1. `using ParallelGraphing`

Functions Currently Supported
1. `make_simple_weighted_graph(size)` - makes a graph of the given size
1. `prims_sequential(graph)` 
1. `prims_parallel(graph)`
1. `dijkstra_all_sources_sequential(graph)`
1. `dijkstra_all_sources_parallel(graph)`
1. `bfs_sequential(graph,source_node)`
1. `bfs_parallel(graph,source_node)`
1. `plot(graph)` - plots graph to allow visual inspection for correctness.
1. `benchmark(graph, algorithm)` - uses the benchmarking suit to run a given algorithm on a given graph where algorithm is an integer (1 for Prim's, 2 for Dijkstra and 3 for BFS)
1. `auto_benchmark(algorithm)` - benchmarks a given algorithm on the provided set of sparsely connected graphs.
1. `auto_connected_benchmark(algorithm)` - benchmarks a given algorithm on the provided set of densely connected graphs.

some simple commands and the order you should probably do them in
1. `add SE751-Assignment` (from pkg manager)
1. `using ParallelGraphing`
1. `g = ParallelGraphing.make_simple_weighted_graph(20)` (or whatever size you want)
1. `ParallelGraphing.plot(g)` displays that graph
1. `ParallelGraphing.benchmark(g)` benchmarks that graph (still under construction)



if you're developing you probably want to use revise before you start making changes....
1. `Julia> CD("location")` (where location is the folder containing SE751-Assignment)
1. `1.0>add Revise`
1. `Julia> using Revise`
1. `1.0>dev SE751-Assignment`
1. `Julia>using ParallelGraphing`



Adding a package
1. navigate inside the SE751-Assignment folder
1. in package manager: `activate .`
1. add package you want to add


Loading a graph
1. `using LightGraphs`
1. `using SimpleWeightedGraphs`
1. `g = LightGraphs.loadgraph("src/graphs/ConnectedGraph10.lg", SWGFormat())` (or other name)

