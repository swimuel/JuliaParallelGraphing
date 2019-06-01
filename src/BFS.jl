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
