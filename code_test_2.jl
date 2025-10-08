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

## test de cette focntion avec un simple graphe 
#file = "G_ex_papier.txt"
file = "gg_05_05_a_1.in"
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
# Fonction pour obtenir les voisins d'un ensemble de sommets
# on a besoin de cette liste car on va faire un problème de séparation et il nous les faut pour les contraintes de type cut-based 
# cela nous servira à garantir que chaque classe induit un sous graphe connecté
#puisque notre ocntarinte fait que si deux sommets u et v sont dans la meme classe i, il existe un chemin connecté entre eux à l'intérieur de cette classe. 
# on va créer un séparateur S: si on peut enlever les sommets de S: u et v seraient dans deux composantes différentes
# si aucun sommet de S est dans la classe i, alors u et v ne sont pas connectés et on a une violation de la contrainte 
# ----------------------------
function boundary_neighbors(G::SimpleGraph, C::Vector{Int})

    #on veut savoir si un sommet est dans le vecteur C
    inC = zeros(Bool, nv(G)) 
    for v in C
        inC[v] = true
    end
    # il faut qu'on stocke les voisins qui sont pas dans C
    N = Int[]
    # ici on boucle sur les sommets de C
    for v in C
        # on boucle sur les voisins de v
        for u in neighbors(G, v)
            if !inC[u] # on vérifie si le voisin n'est pas dans C 
                push!(N, u) # on l'ajoute à la liste des voisins 
            end
        end
    end
    return unique(N) #ici on essaie de vérifier qu'il n'y a pas de doublons 
end
println("==============================")
println("==============================")
println("TEST de la fonction boundary_neighbors")
println("==============================")
println("==============================")
println("Voisins de {1,2} : ", boundary_neighbors(G, [1,2]))
println("Voisins de {3,4} : ", boundary_neighbors(G, [3,4]))
println("Voisins de {5} : ", boundary_neighbors(G, [5]))    
println("Voisins de {1,2,3} : ", boundary_neighbors(G, [1,2,3]))

# =========================================================
# Procédure de séparation 
# Vérifie pour chaque paire (u,v) non adjacente si les deux 
# sont dans la même classe et déconnectés dans le sous-graphe.
# Si oui, ajoute : xu,i + xv,i - sum_{z∈N(u)∪N(v)} x[z,i] ≤ 1
# =========================================================
function separation_on_components!(model::Model, x, G::SimpleGraph, V::Vector{Int}, k::Int; thr::Float64=0.5)
    # on peut compter le nombre de contraintes ajoutées au fur et a mesure du process 
    added = 0
    x_val = value.(x)# on peut récupere la solution courante

    # créer le set des arêtes existantes
    # on peut convertir les arêtes en tuples (min, max) pour éviter les doublons
    # exemple si (u,v) appartient à Eset alors ils sont connectés, sinon ce n'est pas le cas
    Eset = Set( (min(e.src,e.dst), max(e.src,e.dst)) for e in edges(G) )

    #ici on boucle sur les classes 
    for i in 1:k
        # sommets "actifs" de la classe i
        S_i = [v for v in V if x_val[v,i] > thr]
        # si un sommet ou 0 sommet n'est actif, pas de pb de connexité
        if length(S_i) <= 1
            continue
        end


        # graphe induit avec uniquement les sommets de S_i
        G_i, _ = induced_subgraph(G, S_i)
        comps = connected_components(G_i) #on prend les composantes de ces sous graphes 

        # dictionnaire sommet -> composante
        # on essaie de vérifier si deux sommets sont dans la meme composante
        comp_id = Dict{Int,Int}()
        for (cid, comp) in enumerate(comps)
            for v in comp
                comp_id[S_i[v]] = cid
            end
        end

        # pour chaque paire non adjacente dans la même classe mais non connectée
        for u in V, v in V
            # On considère toutes les paires de sommets non adjacente
            if u < v && !( (min(u,v), max(u,v)) in Eset )
                # pour le sommets actifs de la classe i
                if x_val[u,i] > thr && x_val[v,i] > thr
                    # si u et v sont dans des composantes diff alors la contrainte est violée
                    if get(comp_id, u, 0) != get(comp_id, v, 0)
                        #ici on prend tous sles voisins de u et v pour faire le cut based
                        Nuv = union(neighbors(G,u), neighbors(G,v))
                        if !isempty(Nuv)
                            @constraint(model, x[u,i] + x[v,i] - sum(x[z,i] for z in Nuv) <= 1)
                            added += 1  #on a ajouté la contrainte donc on incrémente le compteur
                            println("Ajout contrainte pour u=$u, v=$v, classe=$i")
                        end
                    end
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
file = "gg_05_05_a_1.in"
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
        added = separation_on_components!(model, x, G, V, k; thr=thr)
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
    #println("Solution x : ", value.(x))
end

# Appel :
run_separation!(model, x, G, collect(V), k)








