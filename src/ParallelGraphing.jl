module ParallelGraphing

#Library of sequential and parallel graphs

#dependencies
import LightGraphs
import SimpleWeightedGraphs
import DataStructures
import Plots
import GraphRecipes
import Colors
import BenchmarkTools
import Base.push!

using LightGraphs
using SimpleWeightedGraphs
using Distributed
using DataStructures
using GraphRecipes
using Plots
using Colors
using BenchmarkTools
using Base.Threads: @threads, nthreads, Atomic, atomic_add!, atomic_cas!

include("Dijkstra.jl")
include("BFS.jl")
include("Prims.jl")

export bfs_parallel, bfs_sequential,
	prims_sequential, prims_parallel, prims_priority_queue_sequential,
	dijkstra_all_sources_sequential, dijkstra_all_sources_parallel,
	plot, benchmark, make_simple_weighted_graph


#=helper function to create a SimpleWeightedGraph, graph is randomly made.
	=#
function make_simple_weighted_graph(size)
	g = SimpleWeightedGraph(size)
	for i in 2:size
		add_edge!(g,i-1,i,rand(1:10))

		#generate some random junk edges
		source = rand(1:i)
		destination = rand(1:i)
		if !has_edge(g,source,destination) && source!=destination
			add_edge!(g,source,destination,rand(1:10))
		end
	end
	# plot(g,size)
	return g
end


function plot(graph)
	size = nv(graph)
	Plots.display(graphplot(graph,marker = (:rect),markersize = 1.5,linecolor = :red, names = 1:size))
end

function auto_benchmark(alg)
	for i in 1:20
		g = LightGraphs.loadgraph(string("src/graphs/graph",i,".lg"), SWGFormat())
		println(string("Benchmarking Graph ",i))
		benchmark(g, alg)
		println()
	end
end

function auto_connected_benchmark(alg)
	for i in 21:30
		g = LightGraphs.loadgraph(string("src/graphs/graph",i,".lg"), SWGFormat())
		println(string("Benchmarking Graph ",i))
		benchmark(g, alg)
		println()
	end
end

function benchmark(graph, alg)
	if alg == 1
		para = @benchmark prims_parallel($graph) samples=10 seconds=120 evals=1
		seq = @benchmark prims_sequential($graph) samples=10 seconds=120 evals=1
		display_benchmark_results(para,seq)
	elseif alg == 2
		para = @benchmark dijkstra_all_sources_parallel($graph) samples=10 seconds=120 evals=1
		seq = @benchmark dijkstra_all_sources_sequential($graph) samples=10 seconds=120 evals=1
		display_benchmark_results(para, seq)
	elseif alg == 3
		para = @benchmark bfs_parallel($graph, 1) samples=10 seconds=120 evals=1
		seq = @benchmark bfs_sequential($graph, 1) samples=10 seconds=120 evals=1
		display_benchmark_results(para, seq)
	else
		println("invalid algorithm")
	end
	#dump(result)

end

function display_benchmark_results(para, seq)
	println("average parallel time  : ", match(r"\(([^)]+)\)", string(BenchmarkTools.mean(para))).captures[1])
	println("average sequential time: ", match(r"\(([^)]+)\)", string(BenchmarkTools.mean(seq))).captures[1])
	println("total samples: ", length(para.times))
end
#Make a 20 of graphs as a testbed
#Load using loadgraph("graph{number}",SWGFormat())
#Must be using simple weighed graph

function create_test_graphs()
	for i in 1:20
		g  = make_simple_weighted_graph(i*100)
		savegraph(string("graph",i,".lg"),g)
	end

end

#=helper function that takes a SimpleWeightedGraph and outputs an
	adjacency matrix.
	inputs: @graph: SimpleWeightedGraph type input graph
	outputs: @out_matrix: Adjacency matrix
=#
function graph_to_matrix(graph)
	output_matrix = Matrix(undef,nv(graph),nv(graph))

	for source in vertices(graph)
		for destination in vertices(graph)
			if has_edge(graph,source,destination)
				output_matrix[source,destination] = graph.weights[source,destination]
			else
				output_matrix[source,destination] = 0
			end
		end
	end
	return output_matrix
end

#helper function to construct a graph once the parents have been found for prims, takes in the original matrix
#and an array of the parents of each node.
function construct_graph(parent,matrix)
	source = Int64[]
	destination = Int64[]
	weights = Float64[]
	for i in 2:length(parent)
		Base.push!(source,parent[i])
		Base.push!(destination,i)
		Base.push!(weights,matrix[i,parent[i]])
	end
	return SimpleWeightedGraph(source,destination,weights)
end

end # module
