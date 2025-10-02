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
file = "gg_05_05_a_1.in"  
E, W_vect = readWeightedGraph_paper(file)
n = length(W_vect)
V = 1:n
k = 2 # nombre de partitions

# Dictionnaire de poids
w = Dict(i => W_vect[i] for i in V)
println("Poids des sommets : ", w)
println("Poids total : ", sum(values(w)))

# ----------------------------
# Création du graphe
# ----------------------------
G = Graphs.SimpleGraph(n)
for (u,v) in E
    Graphs.add_edge!(G,u,v)
end


# ----------------------------
# formulation du modèle 
# ----------------------------

model = Model(Gurobi.Optimizer)

#variable définition : ici on a une variable binaire x pour tous les sommets et toutes les classes
@variable(model, x[v in V, i in 1:k], Bin)

#contrainte numéro 1
for i in 1:k-1
    @constraint(model, sum(w[v]*x[v,i] for v in V) <= sum(w[v]*x[v,i+1] for v in V))
end 

#contrainte numéro 2
for v in V
    @constraint(model, sum(x[v,i] for i in 1:k) <=1)
end 

#contrainte 3 indique que chaque classe induit un sous graphe connecté --> donc on ne va pas l'inclure tout de usite car on va utiliser l'algo de séparation 

