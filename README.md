# SE751-Assignment

before initializing Julia run the following command in the shell to create n threads at startup
1. `export JULIA_NUM_THREADS=n`

some simple commands and the order you should probably do them in
1. `add SE751-Assignment` (from pkg manager)
1. `using ParallelGraphing`
1. `g = ParallelGraphing.make_simple_weighted_graph(20)` (or whatever size you want)
1. `ParallelGraphing.plot(g)` displays that graph
1. `ParallelGraphing.benchmark(g)` benchmarks that graph (still under construction



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
1. `g = LightGraphs.loadgraph("src\\graphs\\ConnectedGraph10.lg", SWGFormat())` (or other name)

Functions Currently Supported
1. `prims_sequential(graph)`
1. `prims_parallel(graph)`
1. `dijkstra_all_sources_sequential(graph)`
1. `dijkstra_all_sources_parallel(graph)`
1. `bfs_sequential(graph,source_node)`
1. `bfs_parallel(graph,source_node)`
1. `plot(graph)`
