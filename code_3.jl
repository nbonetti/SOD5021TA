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
#file = "G_ex_papier.txt"
#file = "gg_05_05_a_1.in"
file = "gg_10_10_a_1.in"
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
# amélioration du modèle question 2 avec connectivity inequalities 
# ----------------------------


# ----------------------------
#on utilise plus la fonction boundary_neighbors car on calclule les voisins directement avec union(neighbors(G, u), neighbors(G, v))
# ----------------------------


function separation_connectivity!(model::Model, x, G::SimpleGraph, V::Vector{Int}, k::Int; thr=0.5)
    added = 0
    x_val = value.(x)

    # Créer un set pour les arêtes existantes
    Eset = Set( (min(e.src,e.dst), max(e.src,e.dst)) for e in edges(G) )

    for i in 1:k
        # Sommets actifs dans la classe i
        S_i = [v for v in V if x_val[v,i] > thr]
        if length(S_i) <= 1
            continue
        end

        # Vérifier toutes les paires non adjacentes dans la classe i
        for u in S_i, v in S_i
            if u < v && !((min(u,v), max(u,v)) in Eset)
                # Construire le digraphe Di pour la min-cut
                # Ici on fait une approche simplifiée : on prend voisins communs
                # on en calcule plus les composantes mauis on boucle directement sur toutes les paires de sommets non adjacents dans S_i
                Nuv = union(neighbors(G, u), neighbors(G, v))
                if !isempty(Nuv)
                    # Contrainte de type connectivity lifted
                    @constraint(model, x[u,i] + x[v,i] - sum(x[z,i] for z in Nuv) <= 1)
                    added += 1
                    println("Ajout contrainte connectivity liftée pour u=$u, v=$v, classe=$i")
                end
            end
        end
    end

    return added
end

# ----------------------------
# Programme principal
# ----------------------------
#file = "G_ex_papier.txt"
#file = "gg_05_05_a_1.in"
file = "gg_10_10_a_1.in"
E, W_vect = readWeightedGraph_paper(file)
n = length(W_vect)
V = 1:n
k = 2

w = Dict(i => W_vect[i] for i in V)
println("==============================")
println("==============================")
println("Lecture de l'instance depuis le fichier : ", file)
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
# Formulation du modèle Ck(G,w)
# ----------------------------

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 1)

@variable(model, x[v in V, i in 1:k], Bin)

# (1) Ordre des poids
for i in 1:k-1
    @constraint(model, sum(w[v]*x[v,i] for v in V) <= sum(w[v]*x[v,i+1] for v in V))
end

# (2) Un sommet au plus dans une classe
for v in V
    @constraint(model, sum(x[v,i] for i in 1:k) <= 1)
end

# (3) Connectivité — sera ajoutée par séparation

# Objectif
@objective(model, Max, sum(w[v]*x[v,1] for v in V))

println("==============================")
println("==============================")
println("Formulation du modèle terminée")
println("==============================")
println("==============================")
println("Nombre de variables : ", num_variables(model))
println("Objectif : ", objective_function(model))

# ----------------------------
# Optimisation + séparation
# ----------------------------
# on résoud notre problème avec le modèle "primal" sans les contraintes 3 de l'article
optimize!(model) 

function run_separation!(model, x, G, V, k; max_iter=70, thr=0.5)
    iter = 0
    total_added = 0
    while true
        iter += 1
        println("==============================")
        println("==============================")
        println("Itération de séparation $iter")
        println("==============================")
        println("==============================")
        # on applique la procédure de séparation
        # on peut choisir un seuil thr pour déterminer si un sommet est "dans" une classe
        # si une composante ets violée, on ajoute la contrainte au modèle et on peut augmenter added
        added = separation_connectivity!(model, x, G, V, k; thr=thr)
        total_added += added
        println("Contraintes ajoutées : ", added)
        if added == 0 || iter >= max_iter
            println("Aucune contrainte supplémentaire nécessaire. Fin.")
            break
        end
        # on realnce le modèle avec toutes les contraintes en plus
        optimize!(model)
    end
    optimize!(model)
    println("==============================")
    println("==============================")
    println("Résultat final")
    println("==============================")
    println("==============================")
    println("Statut : ", termination_status(model))
    println("Valeur obj : ", objective_value(model))
    println("Nombre total de contraintes ajoutées : ", total_added)
    #println("Solution x : ", value.(x))
end

# Appel :
run_separation!(model, x, G, collect(V), k)

#pour comparer avec la méthode des composantes
t1 = @elapsed run_separation!(model, x, G, collect(V), k)
println("Temps méthode code 3 : ", t1, " secondes")








