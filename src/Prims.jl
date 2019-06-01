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
