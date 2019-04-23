module ParallelGraphing

#dependencies
import LightGraphs
import SimpleWeightedGraphs

using LightGraphs
using SimpleWeightedGraphs
using Distributed

export prims_sequential, prims_parallel


#=helper function to create a connected SimpleWeightedGraph
	=#
function make_Simple_Weighted_Graph(size)
	g = SimpleWeightedGraph(size)
	for i in 2:size
		add_edge!(g,i-1,i,1)
		#add some junk edges which algorthim has to filter out_matrix
		source = rand(1:i)
		destination = rand(1:i)
		if !has_edge(g,source,destination)
			add_edge!(g,source,destination,10)
		end 

	end
	return g 
end

#=helper function that takes a SimpleWeightedGraph and outputs an
	adjacency matrix. Made parallel with a parallel for loop. 
	inputs: @graph: SimpleWeightedGraph type input graph
	outputs: @out_matrix: Adjacency matrix
=#
function graph_to_matrix(graph)
	output_matrix = Matrix(undef,nv(graph),nv(graph))

	#Made parallel.
	Threads.@threads for source in vertices(graph)
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

		#update the sources and keys array
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
		push!(key,Inf)
		push!(mstSet,false)
		push!(parent,0)
	end

	key[1] = 0
	parent[1] = -1

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

end # module

