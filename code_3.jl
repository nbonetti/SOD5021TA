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
# amélioration du modèle question 2 avec crossing inequalities de la section 3.2 de l'article
# ctet methode est censée améliorer la borne inferieure
# on va s'intéresr aux inégalites numérotées (6) de l'article 
#si un ensemble de sommets S est presque entièrement dans une classe i ( sum sur v x[v,i] est élevée) et dans ce cas S ne peut pas etre entierement affectée à i ou j sans que leur frontière ne soit affectée 



function separation_connectivity!(model::Model, x, G::SimpleGraph, V::Vector{Int}, k::Int; thr=0.5)
    added = 0
    x_val = value.(x)
    # La séparation complète nécessite de trouver toutes les faces du graphe planaire.
    # Pour cet exemple, nous allons simplifier en utilisant un cycle arbitraire
    # comme ensemble de sommets F pour illustrer la logique de l'article.

    # ---------------------------------------------------------------------
    # SIMPLIFICATION pour l'implémentation 
    # ---------------------------------------------------------------------
    # Nous allons prendre tous les sommets comme un "cycle de face" fictif pour l'exemple.
    # Idéalement, F_vertices devrait être le contour d'une face valide.
    F_vertices = collect(V)
    f = length(F_vertices) # Taille de la face |F|
    if f < 4
        return 0 # Nécessite au moins 4 sommets
    end

    # ---------------------------------------------------------------------
    #dans l'article, on nous détaille les matrices L et R, on va donc les construire ici
    # L[j, i] = max {x_val[F(j'), i]} pour j' <= j
    # R[j, i] = max {x_val[F(j'), i]} pour j' >= j
    
    L = zeros(f, k) # Matrice L : L[j, i]
    R = zeros(f, k) # Matrice R : R[j, i]

    #pour la matrice L 
    for i in 1:k
        max_val = -Inf
        for j in 1:f
            max_val = max(max_val, x_val[F_vertices[j], i])
            L[j, i] = max_val
        end
    end

    #pour la matrice R 
    for i in 1:k
        max_val = -Inf
        for j in f:-1:1
            max_val = max(max_val, x_val[F_vertices[j], i])
            R[j, i] = max_val
        end
    end

    #construction de la matrice M est plus compliquée : 
    # pour tous les j apparetnant ) F privé de 1, et pour tout i1 et i2 on a :
    # si j=2 --> M(j,i1,i2)= xtilde[F(1),i1] + xtilde[F(3),i2]
    #sinon max(M(j-1,i1,i2), L(j-1,i1)+xtilde[F(j),i2]
    #notre matrice M est donc en 3 dimensions

    M = Dict{Tuple{Int, Int, Int}, Float64}() 
    
    for i1 in 1:k, i2 in 1:k
        if i1 == i2
            continue
        end

        # Initialisation pour j=2
        v1 = F_vertices[1]
        v2 = F_vertices[2]
        M[2, i1, i2] = x_val[v1, i1] + x_val[v2, i2]
        
        # Calcul par récurrence pour j = 3 à f
        for j in 3:f
            # Cas M(j-1, i1, i2)
            prev_M = M[j - 1, i1, i2]

            # Cas L(j-1, i1) + x_val[F(j), i2]
            v_j = F_vertices[j]
            new_val = L[j - 1, i1] + x_val[v_j, i2]
            
            M[j, i1, i2] = max(prev_M, new_val)
        end
    end

    # puis avec tout cela on peut verifier la violation 
    #on check si M(j-1,i1,i2) + xtilde[F(j),i1] + R(j+1,i2) > 3
    # Le sommet F(j) joue le rôle de t1 dans l'inégalité de croisement xs1 + xs2 + xt1 + xt2 <= 3.
    return added
end

























# ----------------------------
# Programme principal
# ----------------------------
file = "G_ex_papier.txt"
#file = "gg_05_05_a_1.in"
#file = "gg_10_10_a_1.in"
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
        # on relance le modèle avec toutes les contraintes en plus
        optimize!(model)
    end
    if termination_status(model) != OPTIMAL
        optimize!(model)
    end
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



run_separation!(model, x, G, collect(V), k) # <- Premier appel inutile
t1 = @elapsed run_separation!(model, x, G, collect(V), k) # <- Mesure le temps du second appel
println("Temps méthode code 3 : ", t1, " secondes")





