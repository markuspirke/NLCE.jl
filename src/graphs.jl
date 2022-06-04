using CSV, DataFrames, DelimitedFiles

struct Cluster
    name::String
    order::Int64
    sites::Int64 #number of sites
    bonds::Vector{Tuple{Int64,Int64}} #bonds between sites
    C::Float64 #embedding factor
end

"""
Takes filepath as input and returns name of graph.
"""
function graph_name(filename::String)
    name = split(filename, "File")[2]

    name
end
"""
Takes filepath as input and returns name location of graph.
"""
function get_location(filename)
    location = split(filename, "BondFile")[1]

    location
end
"""
Takes filepath as input and returns number of sites of the graph.
"""
function number_sites(filename::String)
    df = CSV.read(filename, DataFrame, header=false, limit=1)
    N = df[!, 1][1]

    N
end
"""
Takes filepath as input and returns number of monomorphisms of the graph.
"""
function monomorphism(filename::String)
    filename_split = split(filename, "-")[3]
    filename_split = split(filename_split, "m")[2]
    N = parse(Int, filename_split)

    N
end

"""
Takes filepath as input and returns number of automorphisms of the graph.
"""
function automorphism(filename::String)
    filename_split = split(filename, "-")[4]
    filename_split = split(filename_split, ".")[1]
    filename_split = split(filename_split, "t")[2]
    N = parse(Int, filename_split)

    N
end

"""
Takes filepath as input and returns bonds of the graph.
"""
function bonds(filename::String)
    M = readdlm(filename, Int, skipstart=2, comments=true, comment_char=';')
    M = M .+ 1
    bs = [(M[i, 2], M[i, 3]) for i in 1:length(M[:, 1])]

    bs
end

"""
Takes filepath as input and return the order in pertubation theory of the graph. 
"""
function get_order(filename::String)
    #number of sites
    N = number_sites(filename)
    #bs vector of tuples of bonds
    bs = bonds(filename)
    #vector which counts how many edges a certain site has 
    nodes = zeros(Int64, N)

    for (i, j) in bs
        nodes[i] += 1
        nodes[j] += 1
    end

    number_bonds = length(bs)
    odd_sites = isodd.(nodes)

    #if graph has already zero or only two odd nodes just return the number of bonds
    if sum(odd_sites) == 0 #|| sum( odd_sites ) == 2 would be needed for longi field term
        return number_bonds
    else #now we need to add bonds as long as we have only zero or two odd nodes 
        for _ in 1:number_bonds
            #take first and last node with odd number of edges and set edge between them
            k = findfirst(odd_sites)
            l = findlast(odd_sites)
            new_bond = (k, l)
            push!(bs, new_bond)
            #add +1 one to the two nodes where the edge was added
            nodes[k] += 1
            nodes[l] += 1
            odd_sites = isodd.(nodes)
            #check again whether the graph is fine now
            if sum(odd_sites) == 0 #|| sum( odd_sites ) == 2
                return length(bs)
            end
        end
    end
end

function no_odd(arr) # check that no odd elements are in arr
	for i in arr
		if i % 2 == 1
			return false
		end
	end
	return true
end

function get_order_matthias(filename)
	file = open(filename, "r")
	num_bonds = 0
	graph = []

	#Read first line which encodes number of sites
	line1 = readline(file)
	num_sites = parse(UInt64, line1)

	degrees = [0 for x=1:num_sites] # Initialize array to save degrees of vertices (i.e. the number of neighbours)

	# Second line is just a dummy parameter
	line2 = readline(file)

	lines = readlines(file)

	# Read the rest of the file line by line. Just split each line and put the second and third number into an array. If semicolon is appended remove it.
	# Surely, you can do this much better.
	for line in lines

		if length(line) != 0 && line[1] != '#'
			sites = split(line, ' ')
			
			src =  parse(UInt64, sites[2])
			if last(sites[3]) == ';'
				tar = parse(UInt64, chop(sites[3]))
			else
				tar =  parse(UInt64, sites[3])
			end	
			
			bond = [src, tar]
			degrees[src+1] += 1
			degrees[tar+1] += 1
			num_bonds += 1
			push!(graph, bond)
		end
	end

	# Will be the minimum number of edges/bonds, which have to be duplicated s.t. no vertex of odd degree remains
	m = copy(num_bonds)

	# Use bit representation of numbers from 0 to num_bonds^(2-1) interpret 1 as bond is duplicated and 0 as bond is not duplicated.
	# Dynamic degree holds the resulting degrees. No_odd checks if a solution has been found.


	for i in 0:(2^num_bonds-1) 
		added_bonds = count_ones(i)
		if added_bonds < m
			dynamic_degrees = copy(degrees)
			bits = last(bitstring(i), num_bonds)	
			if added_bonds < m
				for i in 1:num_bonds
					if bits[i] == '1'
						for site in graph[i]				
							dynamic_degrees[site+1]+=1
						end
					end
				end
			end	
			if no_odd(dynamic_degrees) && added_bonds <= m
				m = copy(added_bonds) 
			end
		end
	end

	return m+num_bonds
end


"""
Input filename and the cluster for that filename and return all indices of Subgraphs
which we need to calculate the reduced energies
"""
function subgraph_embeddings(filename::String, cluster::Cluster)

    location = get_location(filename)
    df = CSV.read(string(location, "Subgraphs", cluster.name), DataFrame, header=["name", "em", "aut"])

    names = df[!, 1]
    names_split = split.(names, "-")
    N = parse.(Int, [name[1] for name in names_split]) .+ 2

    embeddings = [(N[i], df[!, 2][i] / df[!, 3][i]) for i in 1:length(N)]

    embeddings
end

"""
Takes filepath and return a datatype Cluster with all properties of the graph.
"""
function create_cluster(filename::String)
    name = graph_name(filename)
    order = get_order_matthias(filename)
    sites = number_sites(filename)
    bs = bonds(filename)

    C = monomorphism(filename) / automorphism(filename)

    cluster = Cluster(name, order, sites, bs, C)
end