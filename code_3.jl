using JuMP, Gurobi
using Graphs

# ----------------------------
# Fonction pour lire les instances
# ----------------------------
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

# ----------------------------
# Lecture de l'instance
# ----------------------------
## test de cette focntion avec un simple graphe 
file = "G_ex_papier.txt"
#file = "gg_05_05_a_1.in"
E, W_vect = readWeightedGraph_paper(file)
n = length(W_vect)
V = 1:n
k = 2

w = Dict(i => W_vect[i] for i in V)
println("==============================")
println("==============================")
println("TEST de la fonction readWeightedGraph_paper")
println("==============================")
println("==============================")


println("==============================")
println("==============================")
println("Lecture de l'instance depuis le fichier  : ", file)
println("==============================")
println("==============================")
println("Arêtes : ", E)
println("Nombre d'arêtes : ", length(E))
println("Poids des arêtes : ", W_vect)

println("Nombre de sommets : ", n)
println("Poids des sommets : ", w)
println("Poids total : ", sum(values(w)))



G = SimpleGraph(n)
for (u, v) in E
    add_edge!(G, u, v)
end
println("==============================")
println("==============================")
println("Création du graphe à partir de l'instance")
println("==============================")
println("==============================")
println("Graphe créé avec ", nv(G), " sommets et ", ne(G), " arêtes.")
println("Liste des arêtes du graphe : ", edges(G))

# ----------------------------
# Lecture de l'instance
# ----------------------------