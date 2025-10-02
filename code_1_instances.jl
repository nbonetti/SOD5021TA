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

# ----------------------
# Lire l'instance
# ----------------------
file = "gg_15_15_a_1.in"  # Remplacer par ton fichier
E, W_vect = readWeightedGraph_paper(file)
n = length(W_vect)
V = 1:n
k = 3 # nombre de partitions

# Créer dictionnaire de poids
w = Dict(i => W_vect[i] for i in V)

println("Poids des sommets : ", w)
Wtot = sum(values(w))
println("Poids total : ", Wtot)

# ----------------------
# Graphe dirigé avec sources
# ----------------------
sources = [Symbol("s$i") for i in 1:k]
Vd = vcat(V, sources)

Ad = []
for (u,v) in E
    push!(Ad,(u,v)); push!(Ad,(v,u))
end
for s in sources
    for v in V
        push!(Ad,(s,v))
    end
end

# ----------------------
# Modèle JuMP
# ----------------------
model = Model(Gurobi.Optimizer)


set_optimizer_attribute(model, "TimeLimit", 900)  
set_optimizer_attribute(model, "MIPGap", 0.001)   # 0.1 % d'écart toléré



@variable(model, f[Ad] >= 0)
@variable(model, y[Ad], Bin)

@objective(model, Max, sum(f[(sources[1],v)] for v in V if (sources[1],v) in Ad))

for i in 1:k-1
    @constraint(model, sum(f[(sources[i],v)] for v in V if (sources[i],v) in Ad) <=
                        sum(f[(sources[i+1],v)] for v in V if (sources[i+1],v) in Ad))
end

for v in V
    inflow = sum(f[(u,v)] for u in Vd if (u,v) in Ad)
    outflow = sum(f[(v,u)] for u in Vd if (v,u) in Ad)
    @constraint(model, inflow - outflow == w[v])
end

for a in Ad
    @constraint(model, f[a] <= Wtot * y[a])
end

for s in sources
    @constraint(model, sum(y[(s,v)] for v in V if (s,v) in Ad) <= 1)
end

for v in V
    @constraint(model, sum(y[(u,v)] for u in Vd if (u,v) in Ad) <= 1)
end

# ----------------------
# Résolution
# ----------------------
optimize!(model)
println("Status = ", termination_status(model))
println("Valeur optimale = ", objective_value(model))

# ----------------------
# Fonction pour reconstruire partition complète
# ----------------------
function reconstruire_partition(s, V, Ad, y)
    partition = Set{Int}()
    file = []

    # commencer avec les arcs activés de la source
    for v in V
        if value(y[(s,v)]) > 0.5
            push!(partition, v)
            push!(file, v)
        end
    end

    # BFS sur les arcs activés entre sommets
    while !isempty(file)
        u = popfirst!(file)
        for v in V
            if v != u && ( (u,v) in Ad && value(y[(u,v)]) > 0.5 ) && !(v in partition)
                push!(partition, v)
                push!(file, v)
            elseif v != u && ( (v,u) in Ad && value(y[(v,u)]) > 0.5 ) && !(v in partition)
                push!(partition, v)
                push!(file, v)
            end
        end
    end

    return sort(collect(partition))
end

# ----------------------
# Extraire toutes les partitions complètes
# ----------------------
partitions = Dict{Symbol, Vector{Int}}()
for s in sources
    partitions[s] = reconstruire_partition(s, V, Ad, y)
end
println("Partitions complètes : ", partitions)

# ----------------------
# Calcul des poids de chaque partition
# ----------------------
poids_partitions = Dict{Symbol, Int}()
for (s, nodes) in partitions
    poids_partitions[s] = sum(w[v] for v in nodes)
end
println("Poids des partitions : ", poids_partitions)
println("Poids minimum = ", minimum(values(poids_partitions)))
println("Poids maximum = ", maximum(values(poids_partitions)))

# ----------------------
# Tracé du graphe
# ----------------------
g = Graphs.SimpleGraph(length(V))
for (u,v) in E
    Graphs.add_edge!(g, u, v)
end

# Palette de couleurs
palette = [colorant"red", colorant"blue", colorant"green", colorant"orange", colorant"purple"]
colors_index = [0 for _ in V]
for (i,s) in enumerate(sources)
    for v in partitions[s]
        colors_index[v] = i
    end
end

node_colors = [c==0 ? colorant"gray" : palette[c] for c in colors_index]

GraphPlot.gplot(g, nodefillc=node_colors)
