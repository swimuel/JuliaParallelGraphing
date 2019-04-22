module ParallelGraphing

#dependencies
import LightGraphs
import SimpleWeightedGraphs

using LightGraphs
using SimpleWeightedGraphs
using Distributed

export prims_parallel

#=helper function that takes a SimpleWeightedGraph and outputs an
	adjacency matrix
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

function construct_graph(parent,matrix)
	source = Int64[]
	destination = Int64[]
	weights = Float64[]
	for i in 2:length(parent)
		push!(source,parent[i])
		push!(destination,i)
		push!(weights,matrix[i,parent[i]])
	end
	println(source)
	println(destination)
	println(weights)
	return SimpleWeightedGraph(source,destination,weights)
end

function prims_parallel(graph)
	if !is_connected(graph)
		throw("Graph not connected")
	end
	
	matrix = graph_to_matrix(graph)
	#from the matrix, we shall do prims 
	
	parent = []
	key = []
	mstSet = []
	

	#initialize the array
	for i in 1:nv(graph)
		push!(key,Inf)
		push!(mstSet,false)
		push!(parent,0)
	end

	key[1] = 0
	parent[1] = -1

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


end # module

