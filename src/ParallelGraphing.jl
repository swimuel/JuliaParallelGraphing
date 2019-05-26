module ParallelGraphing

#Library of sequential and parallel graphs

#dependencies
import LightGraphs
import SimpleWeightedGraphs
import DataStructures

using LightGraphs
using SimpleWeightedGraphs
using Distributed
using DataStructures

export prims_sequential, prims_parallel, dijkstra_all_sources_sequential, dijkstra_all_sources_parallel, prims_priority_queue_sequential


#=helper function to create a connected SimpleWeightedGraph
	=#
function make_simple_weighted_graph(size)
	g = SimpleWeightedGraph(size)
	for i in 2:size
		add_edge!(g,i-1,i,1)


		#generate some random junk edges
		source = rand(1:i)
		destination = rand(1:i)
		if !has_edge(g,source,destination) && source!=destination
			add_edge!(g,source,destination,rand(2:100))
		end

	end
	return g
end

#Make a dictionary of graphs as a testbed
#Load using loadgraph("graph{number}",SWGFormat())
#Must be using simple weighed graph

function create_test_graphs()
	for i in 1:20
		g  = make_simple_weighted_graph(i*100)
		savegraph(string("graph",i,".lg"),g)
	end
end

function benchmark()
	for i in 1:20
		g = loadgraph(string("graph",i,".lg"),SWGFormat())
		@time prims_sequential(g)
	end
end


#=helper function that takes a SimpleWeightedGraph and outputs an
	adjacency matrix. Made parallel with a parallel for loop.
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


#helper function to construct a graph once the parents have been found for prims
function construct_graph(parent,matrix)
	source = Int64[]
	destination = Int64[]
	weights = Float64[]
	for i in 2:length(parent)
		push!(source,parent[i])
		push!(destination,i)
		push!(weights,matrix[i,parent[i]])
	end
	return SimpleWeightedGraph(source,destination,weights)
end


#=Sequential Variant of the prims algorithm to find the minimum spanning tree
inputs: @graph: SimpleWeightedGraph input, has to be connected
#SLOW VERSION
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
		push!(key,Inf)
		push!(mstSet,false)
		push!(parent,0)
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
	return [Edge{Int64}(parent[v],v) for v in vertices(graph) if parent[v] != 0]
end


#parallelized Variant of the prims algorithms, uses the Threads.@threads construct
#SLOW VERSION

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
		push!(key,Inf)
		push!(mstSet,false)
		push!(parent,0)
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
	return [Edge{Int64}(parent[v],v) for v in vertices(graph) if parent[v] != 0]
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


function prims_priority_queue_sequential(graph)
	num_vertices=nv(graph)
	dist_matrix= weights(graph)

	priority_queue = PriorityQueue{Int64,Float64}() #priority_queue to find min

	finished = zeros(Bool,num_vertices) #has the node been visited
	wt = fill(Inf,num_vertices) #distance matrix
	parents = zeros(Int64,num_vertices)

	priority_queue[1] = 0
	wt[1] = 0

	while !isempty(priority_queue)
		v = dequeue!(priority_queue) #Take the smallest nearest node
		finished[v] = true #this node is now done

		#update step
		for u in neighbors(graph,v)
			finished[u] && continue
			if wt[u] > dist_matrix[u,v]
				wt[u] = dist_matrix[u,v]
				priority_queue[u] = wt[u] #Edges are reranked
				parents[u] = v
			end
		end
	end

	return [Edge{Int64}(parents[v],v) for v in vertices(graph) if parents[v] != 0] #Rerig this to include weights
end

end # module
