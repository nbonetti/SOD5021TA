include("utils.jl")
using JuMP, Gurobi
using Graphs, GraphPlot, Colors

# ----------------------
# Création du modele de flot asymmetrique
# ----------------------
function create_asymmetric_flow_model(E_list, W, Wtot, n, V, k)
    # ----------------------
    # Création du digraphe D
    # ----------------------
    # We only have one source in this asymmetric flow model
    sources = [n+1]

    V_and_sources = vcat(V, sources)
    D = createDigraphFromGraph(E_list, W, sources)
    A = collect(edges(D))
    println("===== Creation du digraph =====")
    println("Arcs du digraphe (A)        : ", A)
    println("Nombre d'arcs |A|           : ", length(A))
    println("Nombre de sommets dans D     : ", nv(D))
    println("Sources (sommets fictifs)   : ", sources)
    println("===============================\n\n")

    # ----------------------
    # Modèle JuMP
    # ----------------------
    model = Model(Gurobi.Optimizer)

    @variable(model, f[e in A, i in 1:k] >= 0)
    @variable(model, y[e in A, i in 1:k], Bin)

    @objective(model, Max, sum(W[v] * sum(y[Edge(u, v), 1] for u in V_and_sources if Edge(u, v) in A) for v in V))

    ## Constraints
    # (14)
    for i in 1:k-1
        @constraint(model, sum(W[v] * sum(y[Edge(u, v), i] for u in V_and_sources if Edge(u, v) in A) for v in V) <=
                            sum(W[v] * sum(y[Edge(u, v), i+1] for u in V_and_sources if Edge(u, v) in A) for v in V))
    end

    # (15)
    for i in 1:k
        @constraint(model, sum(y[Edge(n+1, v), i] for v in V if Edge(n+1, v) in A) <= 1)
    end

    # (16)
    for v in V
        @constraint(model,sum(y[Edge(u, v), i] for i in 1:k, u in V_and_sources if Edge(u, v) in A) <= 1)
    end

    #(17)
    for u in V,v in V
        if u<v #on choisit l'index des sommets pour définir la relation d'ordre
            for i in 1:k
                if Edge(n+1, v) in A
                    @constraint(model, y[Edge(n+1,v),i] + sum(y[Edge(e,u),i] for e in V_and_sources if Edge(e,u) in A) <= 1)
                end 
            end
        end
    end


    #(18)
    for a in A, i in 1:k
        @constraint(model, f[a, i] <= Wtot * y[a, i]) 
    end

    # (19)
    for v in V, i in 1:k
        inflow = sum(f[Edge(u,v), i] for u in V_and_sources if Edge(u,v) in A)
        outflow = sum(f[Edge(v,u), i] for u in V_and_sources if Edge(v,u) in A)
        @constraint(model, outflow - inflow <= 0)
    end

    # (20) 
    for v in V
        total_in = sum(f[Edge(u, v), i] for i in 1:k, u in V_and_sources if Edge(u, v) in A)
        total_out = sum(f[Edge(v, u), i] for i in 1:k, u in V_and_sources if Edge(v, u) in A)
        @constraint(model, total_in - total_out == 1)    
    end

    return model, A, V, y
end


function solve_asymmetric_flow_model(model, y)
    # ----------------------
    # Résolution
    # ----------------------
    JuMP.optimize!(model)
    return model, y
end


function display_asymmetric_flow_results(model, y, A, W, k)
    # Model Solving
    println("Status = ", JuMP.termination_status(model))
    println("Valeur optimale = ", JuMP.objective_value(model))
    println("Temps de résolution (s) = ", JuMP.solve_time(model))


    # Dictionnaire pour stocker le résultat : Sommet -> Classe_Index (V -> i)
    partition = Dict{Int, Int}()
    wp = zeros(Int, k)  # Poids par classe
    minw = 0 

    V = 1:n
    V_and_sources =  1:n+1
    
    # Stocke les variables y* optimisées pour une lecture plus rapide
    y_val = value.(y) 


    # Itération sur chaque sommet v de l'ensemble V (les sommets réguliers)
    for v in V
        is_assigned = false
        
        for i in 1:k
            
            is_class_i = sum(y_val[Edge(u, v), i] for u in V_and_sources if Edge(u, v) in A)
            
            # Vérification de l'appartenance à la classe i
            # L'appartenance est établie si la somme est égale à 1.0 
            # tolérance en cas de solution légèrement non entière
            if isapprox(is_class_i, 1.0, atol=1e-5)
                # Le sommet v appartient à la classe i
                partition[v] = i
                wp[i] += W[v]  
                is_assigned = true
                break 
            end
        end

        minw = minimum(wp)

        if !is_assigned
            @warn "Le sommet $v n'a été assigné à aucune des $k classes. Vérifiez la solution."
            partition[v] = 0 
        end
    end

    classes = [findall(v -> partition[v] == i, 1:n) for i in 1:k]

    println("\n\n===== Résultat final =====")
    println("Nombres de classes souhaitées : ", k)
    println("Partitions : ", classes)
    println("Poids par classe : ", wp, "  -> min = ", minw)
    println("===========================\n\n")
    
    return partition
end


function run_asymmetric_flow_model(E_list, W, Wtot, n, V, k)
    model, A, V, y = create_asymmetric_flow_model(E_list, W, Wtot, n, V, k)
    model, y = solve_flow_model(model, y)
    display_asymmetric_flow_results(model, y, A, W, k)
end






