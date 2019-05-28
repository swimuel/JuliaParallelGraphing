# SE751-Assignment

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
