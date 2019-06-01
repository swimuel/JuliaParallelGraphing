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
