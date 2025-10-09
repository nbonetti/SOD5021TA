# ----------------------
# Fonction pour lire les instances de l'article
# ----------------------
function readWeightedGraph_paper(file::String)
    data = readlines(file)
    
    # Lecture nombre de sommets et d'arêtes
    line = split(data[2], " ")
    n = parse(Int, line[1])
    m = parse(Int, line[2])

    # Initialisation
    W = zeros(Int, n)
    E_list = []

    # Lecture des poids
    for i in 1:n
        line = split(data[3 + i], " ")
        W[i] = parse(Int, line[4])  # 4ème colonne = poids
    end

    # Lecture des arêtes
    for i in 1:m
        line = split(data[4 + n + i], " ")
        orig = parse(Int, line[1])
        dest = parse(Int, line[2])
        push!(E_list, (orig+1, dest+1))  # +1 pour Julia
    end

    return E_list, W
end


# ----------------------
# Fonction pour transformer un graphe non orienté en graphe orienté avec k sources (question1/partie4 de l'article)
# ----------------------
function createDigraphFromGraph(E_list::Vector{Any}, W::Vector{Int64}, sources::Vector{Int64})
    n = length(W)
    k = length(sources)
    D = DiGraph(n + k)
    # Add the arcs on both directions
    for (u,v) in E_list
        add_edge!(D, u, v)
        add_edge!(D, v, u)
    end
    # Add the k sources (all the sources are linked to all the nodes)
    for i in sources
        for v in 1:n
            add_edge!(D, i, v)
        end
    end

    return D
end 