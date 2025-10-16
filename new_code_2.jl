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
# Procédure de séparation (le cœur du problème)
# on veut résoudre par un algorithme de branch and cut la formulation de la section 2
#pour ce faire on va former plusieurs fonctions : 
# - find_violating_cut : qui va chercher une contrainte violée dans une solution fractionnaire
# - separation_callback : qui va appeler find_violating_cut pour chaque paire de sommets non adjacents
# - add_cut_constraints : qui va ajouter les contraintes de coupe au modèle
#- solve_bcp_section2 : qui va créer le modèle, définir les variables, l'objectif, les contraintes de base et lancer le B&C
# ------------------------------------------------------------------

"""
    find_violating_cut(G, x_val, i, u, v)
1.C'est la procédure de séparation qui résout le problème de séparateur de poids minimal (min cut)

2.Pour une classi i, l'inégalité de connexité (3) est violée pour la paire (u,v) si la somme des variables x[z,i], 
pour les sommets z du séparateur minimale est trop petite pour blqouer le chemin entre u et V. 

3. L'article montre que ce problème se ramène à truver la coupe nodale de poids minimum entre u et v dan sun graphe Di.
Le graphe auxiliaire est construit en "dédoublant" chaque sommet `w` en `w_in` et `w_out`.
- Un arc `(w_in -> w_out)` est créé avec une capacité `x_val[w, i]`. C'est la capacité du sommet.
- Pour chaque arête `{a, b}` dans G, on ajoute deux arcs de capacité infinie :
  `(a_out -> b_in)` et `(b_out -> a_in)`.

Le flot max de `u_out` à `v_in` dans ce graphe est égal à la coupe nodale minimale.
"""
function find_violating_cut(G::SimpleGraph, x_val::Dict, i::Int, u::Int, v::Int)
    n = nv(G)
    
    # Indices pour les sommets dédoublés dans le modèle de coupe min
    # w_in -> w, w_out -> n + w
    source_node = n + u # la source est u_out
    sink_node = v     #  le puits est v_in

    flow_model = Model(Gurobi.Optimizer)
    set_silent(flow_model)

    # Variables de flot f[i, j] pour le flot de i à j
    @variable(flow_model, f[1:(2n), 1:(2n)] >= 0)

    # 1. Contraintes de capacité des sommets (arcs w_in -> w_out)
    for w in 1:n
        @constraint(flow_model, f[w, n + w] <= x_val[w, i])
    end

    # 2. Contraintes de capacité des arêtes originales (infinie)
    # C'est une grande valeur, mais suffisante pour ne pas être une contrainte limitante.
    # Les arcs reliant les sommets (a_out -> b_in et b_out -> a_in) ont une capacité "infinie".
    # Cela garantit que la coupe nodale ne se fera JAMAIS sur ces arcs, mais toujours sur les arcs de sommets.
    INF_CAPACITY = sum(w for w in values(x_val)) + 1.0
    for edge in edges(G)
        a, b = src(edge), dst(edge)
        # Arc a_out -> b_in
        @constraint(flow_model, f[n + a, b] <= INF_CAPACITY)
        # Arc b_out -> a_in
        @constraint(flow_model, f[n + b, a] <= INF_CAPACITY)
    end

    # 3. Conservation du flot pour chaque sommet
    # Le flot entrant doit être égal au flot sortant pour tous les nœuds intermédiaires.
    for w in 1:(2n)
        if w != source_node && w != sink_node
            incoming_flow = sum(f[j, w] for j in 1:(2n))
            outgoing_flow = sum(f[w, j] for j in 1:(2n))
            @constraint(flow_model, incoming_flow == outgoing_flow)
        end
    end

    # 4. Objectif : maximiser le flot sortant de la source
    @objective(flow_model, Max, sum(f[source_node, j] for j in 1:(2n)))
    
    optimize!(flow_model)

    # Si le solveur échoue, on ne fait rien
    if termination_status(flow_model) != MOI.OPTIMAL
        return 0.0, Int[]
    end

    min_cut_value = objective_value(flow_model)

    # 5. Vérification de la violation
    # La contrainte est violée si x_u,i + x_v,i - 1 > min_cut_value
    # Si cette condition est vraie, la coupe S trouvée est une contrainte violée (lazy cut)
    if x_val[u, i] + x_val[v, i] - 1.0 > min_cut_value + 1e-6
        # --- Extraction du séparateur S ---
        # Le séparateur S est obtenu en trouvant tous les sommets w du graphe original
        # tels que w_in est atteignable depuis la source dans le graphe résiduel
        # mais w_out ne l'est pas. C'est l'ensemble des sommets dont les arcs (w_in -> w_out)
        # sont saturés par le flot maximum.
        # On trouve les sommets accessibles depuis la source dans le graphe RÉSIDUEL.
        f_val = value.(f)
        
        # Graphe résiduel implicite
        q = [source_node] # File pour le parcours en largeur (BFS)
        reachable = falses(2n)
        reachable[source_node] = true
        
        head = 1
        while head <= length(q)
            curr = q[head]
            head += 1
            
            for next_node in 1:(2n)
                if !reachable[next_node]
                    # Si un flot peut encore passer de curr à next_node
                    # Arc avant : flot < capacité
                    if curr <= n # arc in -> out
                        cap = x_val[curr, i]
                    else # arc out -> in
                        cap = INF_CAPACITY
                    end
                    if f_val[curr, next_node] < cap - 1e-6
                        reachable[next_node] = true
                        push!(q, next_node)
                    end
                    
                    # Arc arrière : flot > 0
                    if f_val[next_node, curr] > 1e-6
                        reachable[next_node] = true
                        push!(q, next_node)
                    end
                end
            end
        end
        
        # Le séparateur S est l'ensemble des sommets `w` du graphe original
        # tels que `w_in` est atteignable mais `w_out` ne l'est pas.
        separator_S = Int[]
        for w in 1:n
            if reachable[w] && !reachable[n + w]
                push!(separator_S, w)
            end
        end
        return min_cut_value, separator_S
    end
    
    # Pas de violation trouvée
    return min_cut_value, Int[]
end


"""
Callback de séparation pour les contraintes de connexité.

La formulation du problème initiale est dite exponentielle car elle contient les inégalités (3)
On ne peut donc pas générer toutes les contraintes à l'avance.
Cette fonction est appelée par le solveur Gurobi à chaque fois qu'une solution fractionnaire
"""
function separation_callback(cb_data, x, G, k)
    n = nv(G)
    V = 1:n
    
    # On ne lance la séparation que sur les solutions fractionnaires aux noeuds du B&C
    # mais Gurobi gère bien avec LazyConstraints.
    
    x_val = Dict((v, i) => callback_value(cb_data, x[v, i]) for v in V, i in 1:k)

    #Itération sur toutes les paires non-adjacentes {u, v} pour chaque classe i :
    for i in 1:k
        for u in V
            for v in (u+1):n
                # On ne teste que les paires non-adjacentes
                if !has_edge(G, u, v)
                    # Heuristique pour ne pas lancer le min-cut pour rien
                    # Si x_u,i et x_v,i sont tous deux proches de 1, cela signifie que u et v
                    # sont susceptibles d'être dans la même classe i, mais qu'il n'y a pas
                    # d'arête entre eux, ce qui pourrait violer la connexité.
                    if x_val[u, i] + x_val[v, i] > 1.0 + 1e-6
                        
                        min_cut_val, S = find_violating_cut(G, x_val, i, u, v)
                        
                        # Si S n'est pas vide, une contrainte violée a été trouvée
                        if !isempty(S)
                            
                            # Ajout de la contrainte (lazy cut)
                            # Le cut est ajouté sous forme de contrainte paresseuse (Lazy Constraint).
                            # Le solveur Gurobi garantit l'optimalité en l'ajoutant à l'arbre B&C.
                            cut = @build_constraint(x[u,i] + x[v,i] - sum(x[z,i] for z in S) <= 1)

                            # MOI.submit : C'est la méthode standard de JuMP/MOI pour 
                            # transmettre le cut au solveur via l'objet de données du callback.
                            MOI.submit(JuMP.backend(model), MOI.LazyConstraint(cb_data), cut)
                            # On peut retourner après avoir ajouté une ou plusieurs coupes
                        end
                    end
                end
            end
        end
    end
end


# ------------------------------------------------------------------
# Programme principal
# ------------------------------------------------------------------
"""
    solve_bcp_section2(instance_file::String, k::Int)
    Cette fonction va construire le problème maitre décrit dans la solve_bcp_section2
"""
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
    
    # --- Symétrie ---
    # Pour casser la symétrie, on peut ordonner les poids des classes
    # C'est l'inégalité (1) de l'article, mais elle est optionnelle si on maximise W_min
    # for i in 1:k-1
    #     @constraint(model, sum(w[v]*x[v,i] for v in V) <= sum(w[v]*x[v,i+1] for v in V))
    # end
    
    println("Modèle initial créé. Lancement du Branch-and-Cut...")


    # Définir une limite de temps globale pour le solveur
    TIME_LIMIT_SECONDS = 100.0
    set_time_limit_sec(model, TIME_LIMIT_SECONDS)


    # set_attribute(model, "LazyConstraints", 1) : Active l'ajout dynamique de contraintes 
    # (Lazy Constraints) pendant l'arbre de Branch-and-Bound.
    set_attribute(model, "LazyConstraints", 1)
    
        start_time = time()

    # set_attribute(model, MOI.LazyConstraintCallback(), ...) : 
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

file ="gg_05_05_a_1.in"
k=2
solve_bcp_section2(file, k)