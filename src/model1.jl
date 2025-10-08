# ======================
# BCPk - Balanced Connected Partition (flow-based)
# Lecture d'instances, partition complète, poids et graphe
# ======================

using JuMP, Gurobi
using Graphs, GraphPlot, Colors

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


# ----------------------
# Lire l'instance
# ----------------------
file = "G_ex_papier.txt"
E_list, W = readWeightedGraph_paper(file)
println("===== Recuperation de l'instance =====")
println("Instance lue ", file)
println("Arêtes du graphe (E)        : ", E_list)
println("Nombre d'arêtes |E|         : ", length(E_list))
println("Poids des sommets (W)       : ", W)
println("===============================\n\n")



# ----------------------
# Initialiser les variables du problème
# ----------------------
Wtot = sum(W)
n = length(W)
V = 1:n
k = 2 
println("===== Problème initialisé =====")
println("Nombre de sommets (n)      : ", n)
println("Nombre de partitions (k)    : ", k)
println("Sommets (V)                 : ", V)
println("Poids des sommets (W)       : ", W)
println("Poids total (Wtot)          : ", Wtot)
println("===============================\n\n")




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

    return model, D, A, sources, V, V_and_sources
end



# # ----------------------
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
