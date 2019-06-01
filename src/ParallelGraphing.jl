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

export bfs_parallel, bfs_sequential ,prims_sequential, prims_parallel, dijkstra_all_sources_sequential, dijkstra_all_sources_parallel, prims_priority_queue_sequential, plot


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

function benchmark(graph, x)
	if x == 1
		println("running prims in parallel")
		display_benchmark_results(@benchmark prims_parallel($graph) samples=10 seconds=120 evals=1)
	elseif x == 2
		println("running prims sequentially")
		display_benchmark_results(@benchmark prims_sequential($graph) samples=10 seconds=120 evals=1)
	else
		println("invalid algorithm")
	end
	#dump(result)

end

function display_benchmark_results(result)
	println("min: ", minimum(result))
	println("median: ", BenchmarkTools.median(result))
	println("mean: ", BenchmarkTools.mean(result))
	println("max: ", maximum(result))
	println("total time: ", sum(result.times) / 1000, "μs")
	println("total samples: ", length(result.times))
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


#threadsafe queue data storage
struct ThreadQueue{T,N<:Integer}
    data::Vector{T}
    head::Threads.Atomic{N} #Index of the head
    tail::Threads.Atomic{N} #Index of the tail
end



#helper frunctions for managing thread queue
function ThreadQueue(T::Type, maxlength::N) where N <: Integer
    q = ThreadQueue(Vector{T}(undef, maxlength), Atomic{N}(1), Atomic{N}(1))
    return q
end

function push!(q::ThreadQueue{T,N}, val::T) where T where N
    # TODO: check that head > tail
    offset = atomic_add!(q.tail, one(N))
    q.data[offset] = val
    return offset
end

function popfirst!(q::ThreadQueue{T,N}) where T where N
    # TODO: check that head < tail
    offset = atomic_add!(q.head, one(N))
    return q.data[offset]
end

function isempty(q::ThreadQueue{T,N}) where T where N
    return (q.head[] == q.tail[]) && q.head != one(N)
end

#iterates through the queue and adds newly found nodes to queue
function bfs_parallel_main(
        next::ThreadQueue,
        g::SimpleWeightedGraph,
        parents::Array{Threads.Atomic{T}},
        level::Array{T}
    ) where T <: Integer
    Threads.@threads for src in level
        vertexneighbors = neighbors(g, src)
        for vertex in vertexneighbors
            parent = atomic_cas!(parents[vertex], zero(T), src)
            if parent == 0
                push!(next, vertex)
            end
        end
    end
end

#initalizes the array of parents and calls the main iteration loop
function bfs_parallel!(
        next::ThreadQueue,
        g::SimpleWeightedGraph,
        source::T,
        parents::Array{Threads.Atomic{T}}
    ) where T<:Integer
    parents[source][] = source
    push!(next, source)
    while !isempty(next)
        level = next.data[next.head[]:(next.tail[] - 1)]
        next.head[] = next.tail[]
        bfs_parallel_main(next, g, parents, level)
    end
    return parents
end

#creates thread queue to iterate through and calls the bfs_parallel! function
function bfs_parallel(g::SimpleWeightedGraph, source::T, nv::T) where T <: Integer
    next = ThreadQueue(T, nv)
    parents = [Atomic{T}(0) for i = 1:nv]
    ParallelGraphing.bfs_parallel!(next, g, source, parents)
	i = [i[] for i in parents]
    LightGraphs.tree([i[] for i in parents]) #create tree from parent node array
end
#finds the number of verticies in the graph then calls the above function using that value
function bfs_parallel(g::SimpleWeightedGraph, source::T) where T <: Integer
    nvg = nv(g)
    ParallelGraphing.bfs_parallel(g, source, nvg)
end

#iterates through all the neighbors assigning them values and then looking for their neighbors
function bfs_sequential_main(next::Array, g::SimpleWeightedGraph, parents::Array, start::Integer)
	for i in next
		if  parents[i] == 0
			parents[i] = start
		end
	end
	for i in next
		if  parents[i] == start
			forward = neighbors(g,i)
			if !(Base.isempty(forward))
				bfs_sequential_main(forward, g, parents, i)
			end
		end
	end
end

#creates next array to iterate through and calls the bfs_sequential_main function
function bfs_sequential(g::SimpleWeightedGraph, source::T, nv::T) where T <: Integer
    next = neighbors(g, source)
    parents = [0 for i = 1:nv]
	parents[source] = source
	#while (0 in parents)
	bfs_sequential_main(next, g, parents, source)
	#end
	LightGraphs.tree([i for i in parents])  #create tree from parent node array
end

#finds the number of verticies in the graph then calls the above function using that value
function bfs_sequential(g::SimpleWeightedGraph, source::T) where T <: Integer
	nvg = nv(g)
	ParallelGraphing.bfs_sequential(g, source, nvg)
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


#=Sequential Variant of the prims algorithm to find the minimum spanning tree
inputs: @graph: SimpleWeightedGraph input, has to be connected

=#
function prims_sequential(graph)
	if !is_connected(graph)
		throw("Graph not connected")
	end

	matrix = graph_to_matrix(graph)
	#from the matrix, we shall do prims

	parent = [] #Sources array
	key = [] #Weights array
	mstSet = [] #Whether or not element is in the minimum spanning tree

	#initialize the array
	for i in 1:nv(graph)
		Base.push!(key,Inf)
		Base.push!(mstSet,false)
		Base.push!(parent,0)
	end

	key[1] = 0
	parent[1] = 0

	for vertex in 1:nv(graph)
		#pick minimum key vertex from set of vertices not in mst
		min = Inf
		min_index = -1

		for v in 1:nv(graph)
			if mstSet[v] ==false && key[v] < min
				min = key[v]
				min_index = v
			end
		end

		#add picked vertex to the MST
		mstSet[min_index] = true

		#update the sources and keys array
		for i in 1:nv(graph)
			if matrix[min_index,i]!=0 && mstSet[i] == false && matrix[min_index,i] < key[i]
				parent[i] = min_index
				key[i] = matrix[min_index,i]
			end
		end
	end
	#format the output
	graph = construct_graph(parent,matrix)
	return graph
end


#parallelized Variant of the prims algorithms, uses the Threads.@threads construct

function prims_parallel(graph)
	if !is_connected(graph)
		throw("Graph not connected")
	end

	matrix = graph_to_matrix(graph)
	#from the matrix, we shall do prims

	parent = [] #Sources array
	key = [] #Weights array
	mstSet = [] #Whether or not element is in the minimum spanning tree

	#initialize the array
	for i in 1:nv(graph)
		Base.push!(key,Inf)
		Base.push!(mstSet,false)
		Base.push!(parent,0)
	end

	key[1] = 0
	parent[1] = 0

	for vertex in 1:nv(graph)
		#pick minimum key vertex from set of vertices not in mst
		min_reduction_array = zeros(Threads.nthreads())
		fill!(min_reduction_array,Inf)

		index_array = zeros(Threads.nthreads())
		fill!(index_array,-1)

		Threads.@threads for v in 1:nv(graph)
			if mstSet[v] ==false && key[v] < min_reduction_array[Threads.threadid()]
				min_reduction_array[Threads.threadid()] = key[v] #Each thread will write it's max and then it will be reduced manually
				index_array[Threads.threadid()] = v
			end
		end
		min_index = trunc(Int,index_array[argmin(min_reduction_array)])

		#add picked vertex to the MST
		mstSet[min_index] = true

		#update the sources and keys array. can be made parallel but is somehow slower.
		for i in 1:nv(graph)
			if matrix[min_index,i]!=0 && mstSet[i] == false && matrix[min_index,i] < key[i]
				parent[i] = min_index
				key[i] = matrix[min_index,i]
			end
		end

	end
	graph = construct_graph(parent,matrix)
	return graph
end


#=dijkstra for a single node
	helper method for all_sources
	takes an adjancency matrix and the source with which is being searched from and returns an array of shortest distances
	=#
function dijk_single(matrix, src)

	dist = Array{Float64}(undef,size(matrix,1)) #keep track of distances to the src
	shortest_path = falses(size(matrix,1)) #Array to keep track of whether a node is in the shortest path

	fill!(dist,Inf)
	fill!(shortest_path,false)

	dist[src] = 0
	for i in 1:size(matrix,1)

		min = Inf
		min_index = -1

		#find minimum distance node
		for i in 1:size(matrix,1)
			if !shortest_path[i] && dist[i] < min
				min = dist[i]
				min_index = i
			end
		end

		shortest_path[min_index] = true

		#update distances
		for i in 1:size(matrix,1)
			if !shortest_path[i] && matrix[min_index,i] != 0 && dist[min_index] != Inf && dist[i] > dist[min_index] + matrix[min_index,i]
				dist[i] =dist[min_index] + matrix[min_index,i]
			end
		end
	end
	return dist
end

#=dijkstra algorthim to find the shortest path from all nodes to all nodes
inputs: @graph SimmpleWeightedGraph, has to be connected
output: Array containing arrays of shortest distances for each node in the graph
=#
function dijkstra_all_sources_sequential(graph)
	if !is_connected(graph)
		throw("Graph not connected")
	end
	matrix = graph_to_matrix(graph)
	#from the matrix, we shall do dijkstra

	weight_matrix = Array{Any}(undef,size(matrix,1)) #matrix that will contain all shortest_paths

	for src in 1:size(matrix,1)
		weight_matrix[src] = dijk_single(matrix,src)
	end
	return weight_matrix
end

#=dijkstra algorthim to find the shortest path from all nodes to all nodes
inputs: @graph SimmpleWeightedGraph, has to be connected
output: Array containing arrays of shortest distances for each node in the graph

Uses threads and the number of threads has to be set before julia starts up
=#
function dijkstra_all_sources_parallel(graph)
	if !is_connected(graph)
		throw("Graph not connected")
	end
	matrix = graph_to_matrix(graph)
	#from the matrix, we shall do dijkstra algorthim

	weight_matrix = Array{Any}(undef,size(matrix,1)) #matrix that will contain all shortest_paths

	Threads.@threads for src in 1:size(matrix,1)
		weight_matrix[src] = dijk_single(matrix,src)
	end
	return weight_matrix
end



end # module
