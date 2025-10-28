using JuMP, Gurobi, Graphs
import MathOptInterface as MOI

# ------------------------------------------------------------------
# Fonction de lecture du graphe (inchangée)
# ------------------------------------------------------------------
function readWeightedGraph_paper(file::String)
    data = readlines(file)
    line = split(data[2], " ")
    n = parse(Int, line[1])
    m = parse(Int, line[2])
    W = zeros(Int, n)
    E_list = []
    for i in 1:n
        line = split(data[3 + i], " ")
        W[i] = parse(Int, line[4])
    end
    for i in 1:m
        line = split(data[4 + n + i], " ")
        orig = parse(Int, line[1])
        dest = parse(Int, line[2])
        push!(E_list, (orig + 1, dest + 1))
    end
    g = SimpleGraph(n)
    for (u, v) in E_list
        add_edge!(g, u, v)
    end
    w_dict = Dict(i => W[i] for i in 1:n)
    return g, w_dict
end


# ------------------------------------------------------------------
# amélioration du modèle de sépartaion résolu en question 2
# ------------------------------------------------------------------

"""
Question 2 récap : 

Notre  programme d'avant implémente correctement le modèle de flot maximum/coupe minimum pour trouver le séparateur,
mais il le fait sur le graphe dédoublé complet, sans les optimisations mentionnées dans la Section 3.1 de l'article.

Améliorations proposées : Question 3
La contraction d'arcs est une optimisation qui vise à réduire la taille du modèle de flot (et donc son temps de résolution) 
en exploitant les variables qui ont déjà atteint une valeur entière (0 ou 1) dans la solution fractionnaire courante
1. Le modèle de flot est construit sur un sous-graphe de taille réduite (moins de 2n sommets).
2. Si  xtilde[v,i] =1 le sommet v et ses voisins assurément connectés à lui sont contractés.
3. Si  xtilde[v,i] =0 le sommet v et ses arcs adjacents à la classe i sont exclus du modèle de flot pour la classe i
"""
# -----------------------------------------------------------------
# -----------------------------------------------------------------
#          Fonction find_violating_cut_optimized
# -----------------------------------------------------------------
# -----------------------------------------------------------------


"""
    find_violating_cut_optimized(G, x_val, i, u, v)

Procédure de séparation optimisée par 'Contraction d'arcs' (réduction du graphe de flot).
Le modèle de flot n'est construit que sur les sommets qui sont :
1. Le séparateur (u et v).
2. Les sommets dont la variable x[w, i] est fractionnaire (0 < x < 1).

ceci réduit significativement la taille du PL de flot à chaque itération.
"""


function find_violating_cut_optimized(G::SimpleGraph, x_val::Dict, i::Int, u::Int, v::Int)
    n = nv(G)
    V = 1:n
    
    # 1. IDENTIFICATION DES SOMMETS PERTINENTS (V')
    # Les sommets avec x=0 ou x=1 sont contractés/exclus.
    V_relevant = Int[] # Contient les indices des sommets originaux (sommets retenus)
    node_map = Dict{Int, Int}() # Mappe l'ancien indice w vers le nouvel indice w'(réduit)
    
    # Les sommets avec x=0 ou x=1 sont contractés (c'est-à-dire exclus de la coupe)
    # On garde seulement les sommets fractionnaires + la paire {u, v}
    for w in V
        is_fractional = (0.0 + 1e-6 < x_val[w, i] < 1.0 - 1e-6)
        
        if is_fractional || w == u || w == v
            push!(V_relevant, w)
            # Le nouvel indice sera l'ordre d'apparition dans V_relevant
            node_map[w] = length(V_relevant)
        end
    end
    
    # n_prime est le nombre de sommets dans le graphe de flot réduit (taille V')
    n_prime = length(V_relevant)
    
    # Si le nombre de sommets est trop petit (e.g. u et v ne sont pas connectés via un noeud fractionnaire)
    if n_prime < 2 || node_map[u] == node_map[v]
         return 0.0, Int[]
    end
    

    # Définition de la source et du puits dans le graphe réduit (taille 2 * n_prime)
    # Indices du graphe de flot (nouveaux) : w'_in = w', w'_out = n' + w'
    source_prime = node_map[u] + n_prime # u'_out
    sink_prime = node_map[v]      # v'_in



    flow_model = Model(Gurobi.Optimizer)
    set_silent(flow_model)



    # Le modèle de flot est de taille (2 * n_prime) x (2 * n_prime) : BEAUCOUP PLUS PETIT
    @variable(flow_model, f[1:(2n_prime), 1:(2n_prime)] >= 0)



    # Grande capacité (infinie) pour les arcs d'arêtes
    INF_CAPACITY = sum(w for w in values(x_val)) + 1.0



    # 2. CONSTRUCTION DU MODÈLE DE FLOT RÉDUIT
    
    # 2.1. Contraintes de capacité des sommets (arcs w'_in -> w'_out)
    for w_prime in 1:n_prime
        w = V_relevant[w_prime] # Indice du sommet original
        
        w_in_prime = w_prime
        w_out_prime = w_prime + n_prime
        
        # Capacité = x_val[w, i] (sera 1.0 pour u/v s'ils étaient entiers, et fractionnaire sinon)
        @constraint(flow_model, f[w_in_prime, w_out_prime] <= x_val[w, i])
    end

    # 2.2. Contraintes de capacité des arêtes originales (infinie)
    # On n'ajoute que les arêtes entre les sommets qui sont restés dans V'.
    for w in V_relevant
        w_out_prime = node_map[w] + n_prime
        
        for z in V_relevant
            z_in_prime = node_map[z]
            
            # On vérifie seulement les arêtes entre les sommets inclus dans le modèle réduit
            if has_edge(G, w, z) 
                # Arc w'_out -> z'_in
                @constraint(flow_model, f[w_out_prime, z_in_prime] <= INF_CAPACITY)
            end
        end
    end

    # 2.3. Conservation du flot
    for p in 1:(2n_prime) # p est l'indice du nœud dans le modèle réduit
        if p != source_prime && p != sink_prime
            incoming_flow = sum(f[j, p] for j in 1:(2n_prime))
            outgoing_flow = sum(f[p, j] for j in 1:(2n_prime))
            @constraint(flow_model, incoming_flow == outgoing_flow)
        end
    end

    # 2.4. Objectif : maximiser le flot sortant de la source
    @objective(flow_model, Max, sum(f[source_prime, j] for j in 1:(2n_prime)))
    
    optimize!(flow_model)

    if termination_status(flow_model) != MOI.OPTIMAL
        return 0.0, Int[]
    end

    min_cut_value = objective_value(flow_model)

    # 3. VÉRIFICATION DE LA VIOLATION ET EXTRACTION DU SÉPARATEUR S
    if x_val[u, i] + x_val[v, i] - 1.0 > min_cut_value + 1e-6
        # --- Extraction du séparateur S par BFS sur le graphe résiduel ---
        f_val = value.(f)
        
        q = [source_prime] 
        reachable = falses(2n_prime)
        reachable[source_prime] = true
        
        head = 1
        while head <= length(q)
            curr_prime = q[head]
            head += 1
            
            # 3.1. Trouver tous les voisins (next_node_prime)
            for next_node_prime in 1:(2n_prime)
                if !reachable[next_node_prime]
                    
                    # 3.2. Calcul de la capacité résiduelle (cap - f_val)
                    cap = 0.0
                    
                    # Case A: Arc V_in -> V_out (Capacité de sommet)
                    if curr_prime <= n_prime && next_node_prime == curr_prime + n_prime
                        w = V_relevant[curr_prime]
                        cap = x_val[w, i]
                    # Case B: Arc V_out -> V_in (Capacité d'arête)
                    elseif curr_prime > n_prime && next_node_prime <= n_prime
                        cap = INF_CAPACITY
                    end
                    
                    # Arc avant (Forward arc) : flot < capacité
                    if cap > 0.0 && f_val[curr_prime, next_node_prime] < cap - 1e-6
                        reachable[next_node_prime] = true
                        push!(q, next_node_prime)
                    end
                    
                    # Arc arrière (Backward arc) : flot > 0
                    if f_val[next_node_prime, curr_prime] > 1e-6
                        reachable[next_node_prime] = true
                        push!(q, next_node_prime)
                    end
                end
            end
        end
        
        # 3.3. Extraction de S (dans les indices originaux)
        separator_S = Int[]
        for w_prime in 1:n_prime
            w = V_relevant[w_prime] # On récupère l'indice original
            # Le sommet w est dans la coupe si w_in est atteignable et w_out ne l'est pas.
            if reachable[w_prime] && !reachable[w_prime + n_prime]
                push!(separator_S, w)
            end
        end
        
        return min_cut_value, separator_S
    end
    
    return min_cut_value, Int[]
end



# ------------------------------------------------------------------
# Callback de séparation
# ------------------------------------------------------------------    

function separation_callback(cb_data, x, G, k)
    n = nv(G)
    V = 1:n
    
    # On ne lance la séparation que sur les solutions fractionnaires aux noeuds du B&C
    # mais Gurobi gère bien avec LazyConstraints.
    
    x_val = Dict((v, i) => callback_value(cb_data, x[v, i]) for v in V, i in 1:k)
    
    for i in 1:k
        for u in V
            for v in (u+1):n
                if !has_edge(G, u, v)
                    if x_val[u, i] + x_val[v, i] > 1.0 + 1e-6
                        
                        # CHANGEMENT ICI : Appel à la fonction optimisée
                        min_cut_val, S = find_violating_cut_optimized(G, x_val, i, u, v)
                        
                        if !isempty(S)
                            model = owner_model(x[1, 1])
                            cut = @build_constraint(x[u,i] + x[v,i] - sum(x[z,i] for z in S) <= 1)
                            MOI.submit(JuMP.backend(model), MOI.LazyConstraint(cb_data), cut)
                            return
                        end
                    end
                end
            end
        end
    end
end


# ------------------------------------------------------------------
# Fonction principale pour résoudre le BCP avec l'amélioration
# ------------------------------------------------------------------

function solve_bcp_section2(instance_file::String, k::Int)
    G, w = readWeightedGraph_paper(instance_file)
    n = nv(G)
    V = 1:n
    
    println("==============================")
    println("Instance: $instance_file, k=$k")
    println("Graphe créé avec $n sommets et $(ne(G)) arêtes.")
    println("==============================")
    
    model = Model(Gurobi.Optimizer)
    
    
    @variable(model, x[v in V, i in 1:k], Bin)
    @variable(model, W_min) # Poids de la classe la plus légère

    # L'objectif est de maximiser le poids de la partition de poids minimum
    @objective(model, Max, W_min)
    
    # Contrainte de définition de W_min
    for i in 1:k
        @constraint(model, W_min <= sum(w[v] * x[v, i] for v in V))
    end
    
    # Contrainte de partition : chaque sommet est dans au plus une classe
    # (permet les k-subpartitions, comme suggéré dans l'article)
    for v in V
        @constraint(model, sum(x[v, i] for i in 1:k) <= 1)
    end
    
    
    
    println("Modèle initial créé. Lancement du Branch-and-Cut...")


    # Définir une limite de temps globale pour le solveur
    TIME_LIMIT_SECONDS = 100.0
    set_time_limit_sec(model, TIME_LIMIT_SECONDS)


    # set_attribute(model, "LazyConstraints", 1) : Active l'ajout dynamique de contraintes 
    # (Lazy Constraints) pendant l'arbre de Branch-and-Bound.
    set_attribute(model, "LazyConstraints", 1)
    
        start_time = time()


    # Enregistre la fonction 'separation_callback' comme la routine à exécuter
    # chaque fois que Gurobi trouve une solution entière (ou presque entière)
    # dans le Branch-and-Bound.
    set_attribute(model, MOI.LazyConstraintCallback(), 
        cb_data -> separation_callback(cb_data, x, G, k)
    )

    println("Modèle prêt. Démarrage du Branch-and-Cut avec une limite de temps globale de $(TIME_LIMIT_SECONDS)s...")
    

    
    optimize!(model)
    end_time = time()

    println("==============================")
    println("Résultat Final (Formulation Section 2)")
    println("==============================")
    
    status = termination_status(model)
    println("Statut : ", status)

    if status == MOI.OPTIMAL
        println("SOLUTION OPTIMALE TROUVÉE.")
    elseif status == MOI.TIME_LIMIT
        println("LIMITE DE TEMPS ATTEINTE.")
    else
        println("Le solveur s'est arrêté pour une autre raison.")
    end

    if has_values(model)
        println("Meilleure solution trouvée (poids min) : ", objective_value(model))
        println("Meilleure borne duale : ", objective_bound(model))
        gap = relative_gap(model)
        println("Gap d'optimalité relatif : ", round(100*gap, digits=2), "%")
        
        println("\nDétail de la meilleure partition trouvée :")
        for i in 1:k
            partition_i = [v for v in V if value(x[v,i]) > 0.5]
            weight_i = sum(w[v] for v in partition_i)
            println("Classe $i (poids $weight_i): $partition_i")
        end
    else
        println("Aucune solution réalisable n'a été trouvée.")
    end
    
    println("\nTemps total d'exécution : ", round(end_time - start_time, digits=2), "s")
    println("Nombre de nœuds B&C explorés : ", node_count(model))
end



#type du fichier


#file1 ="gg_05_05_a_1.in"
#file2 ="gg_05_05_a_2.in"
#file3 ="gg_05_05_a_3.in"
#file4 ="gg_05_05_a_4.in"
#file5 ="gg_15_15_a_1.in"
#file6 ="rnd_20_30_a_1.in"
#file7 ="rnd_20_30_a_2.in"
#file8 ="rnd_20_50_a_1.in"


#file9 ="gg_05_05_b_1.in"
#file10 ="rnd_30_50_a_1.in"
#file11 ="rnd_30_50_a_10.in"
#file12 ="rnd_30_50_c_1.in"


#file13 ="rnd_30_70_a_1.in"
#file14 ="rnd_30_70_a_2.in"
file15 ="rnd_30_70_a_3.in"


#nombre de classes
k=2

# appel de la fonction finale 

solve_bcp_section2(file15, k)
