# ======================
# BCPk - Balanced Connected Partition (flow-based)
# Lecture d'instances, partition complète, poids et graphe
# ======================
include("utils.jl")
using JuMP, Gurobi
using Graphs, GraphPlot, Colors

# ----------------------
# Création du modele de flot
# ----------------------
function create_flow_model(E_list, W, Wtot, n, V, k)
    # ----------------------
    # Création du digraphe D
    # ----------------------
    sources = collect(n+1:n+k)
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

    @variable(model, f[e in A] >= 0)
    @variable(model, y[e in A], Bin)

    @objective(model, Max, sum(f[Edge(sources[1],v)] for v in V if Edge(sources[1],v) in A))

    ## Constraints
    # (7)
    for i in 1:k-1
        @constraint(model, sum(f[Edge(sources[i],v)] for v in V if Edge(sources[i],v) in A) <=
                            sum(f[Edge(sources[i+1],v)] for v in V if Edge(sources[i+1],v) in A))
    end

    # (8)
    for v in V
        inflow = sum(f[Edge(u,v)] for u in V_and_sources if Edge(u,v) in A)
        outflow = sum(f[Edge(v,u)] for u in V_and_sources if Edge(v,u) in A)
        @constraint(model, inflow - outflow == W[v])
    end

    # (9)
    for e in A
        @constraint(model, f[e] <= Wtot * y[e]) 
    end

    # (10)
    for s in sources
        @constraint(model, sum(y[Edge(s,v)] for v in V if Edge(s,v) in A) <= 1)
    end

    # # (11)
    for v in V
        @constraint(model, sum(y[Edge(u,v)] for u in  V_and_sources if Edge(u,v) in A) <= 1)
    end

    return model, A, sources, V, y
end


function solve_flow_model(model, y)
    # ----------------------
    # Résolution
    # ----------------------
    JuMP.optimize!(model)
    return model, y
end

function display_results(model, y, A, sources, W, k)
    # Model Solving
    println("Status = ", JuMP.termination_status(model))
    println("Valeur optimale = ", JuMP.objective_value(model))
    println("Temps de résolution (s) = ", JuMP.solve_time(model))

    # # Visualisation des variables
    yval = value.(y)

    # Forêt dirigée des arcs sélectionnés (y = 1)
    H = DiGraph(n + k)
    for e in A
        if yval[e] > 0.5
            add_edge!(H, src(e), dst(e))
        end
    end

    # Affectation par simple exploration depuis chaque source
    part = fill(0, n)                
    for (i, s) in enumerate(sources)
        stack = [s]
        seen = falses(n + k); seen[s] = true
        while !isempty(stack)
            u = pop!(stack)
            for v in outneighbors(H, u)
                if !seen[v]
                    seen[v] = true
                    push!(stack, v)
                    v <= n && (part[v] = i)   # affecter si c'est un vrai sommet
                end
            end
        end
    end

    # Résumé final des partitions trouvées
    classes = [findall(v -> part[v] == i, 1:n) for i in 1:k]
    wp = [sum(W[v] for v in classes[i]) for i in 1:k]
    minw = minimum(wp)

    println("\n\n===== Résultat final =====")
    println("Nombres de classes souhaitées : ", k)
    println("Partitions : ", classes)
    println("Poids par classe : ", wp, "  -> min = ", minw)
    println("===========================\n\n")
end 

function run_flow_model(E_list, W, Wtot, n, V, k)
    model, A, sources, V, y = create_flow_model(E_list, W, Wtot, n, V, k)
    model, y = solve_flow_model(model, y)
    display_results(model, y, A, sources, W, k)
end






