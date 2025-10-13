using JuMP, Gurobi, Graphs
import MathOptInterface as MOI

# ----------------------------
# Fonctions d'aide (inchangées)
# ----------------------------

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
        push!(E_list, (orig+1, dest+1))
    end
    return E_list, W
end

function boundary_neighbors(G::SimpleGraph, C::Vector{Int})
    inC = zeros(Bool, nv(G))
    for v in C
        inC[v] = true
    end
    N = Int[]
    for v in C
        for u in neighbors(G, v)
            if !inC[u]
                push!(N, u)
            end
        end
    end
    return unique(N)
end

# =========================================================
# Procédure de séparation avec CALLBACK pour Branch-and-Cut
# =========================================================
function component_separation_callback(cb_data, x, G, V, k; thr::Float64=0.5)
    added_cuts = 0
    
    try
        x_val = Dict((v, i) => JuMP.callback_value(cb_data, x[v, i]) for v in V, i in 1:k)
        
        for i in 1:k
            S_i = [v for v in V if x_val[v,i] > thr]

            if length(S_i) <= 1
                continue
            end

            G_i, _ = induced_subgraph(G, S_i)
            comps = connected_components(G_i) 

            if length(comps) > 1
                println("Classe $i : ", length(comps), " composantes détectées.")
                for j in 2:length(comps)
                    component = comps[j]
                    C = [S_i[v] for v in component]
                    N_C = boundary_neighbors(G, C)

                    lhs = sum(x[v,i] for v in C) - sum(x[z,i] for z in N_C)
                    rhs = length(C) - 1
                    
                    MOI.submit(model, MOI.UserCut(cb_data), 
                        JuMP.build_constraint(model, lhs <= rhs))
                    
                    added_cuts += 1
                end
            end
        end
    catch
        return 0
    end
    return added_cuts
end

# ----------------------------
# Programme principal
# ----------------------------
file = "gg_10_10_a_1.in"
E, W_vect = readWeightedGraph_paper(file)
n = length(W_vect)
V = 1:n
k = 2
w = Dict(i => W_vect[i] for i in V)
G = SimpleGraph(n)
for (u, v) in E
    add_edge!(G, u, v)
end
println("==============================")
println("Création du graphe à partir de l'instance")
println("==============================")
println("Graphe créé avec ", nv(G), " sommets et ", ne(G), " arêtes.")


# ----------------------------
# Formulation du modèle Ck(G,w)
# ----------------------------

model = Model(Gurobi.Optimizer)
@variable(model, x[v in V, i in 1:k], Bin)

# (1) Ordre des poids
for i in 1:k-1
    @constraint(model, sum(w[v]*x[v,i] for v in V) <= sum(w[v]*x[v,i+1] for v in V))
end

# (2) Un sommet au plus dans une classe
for v in V
    @constraint(model, sum(x[v,i] for i in 1:k) <= 1)
end

# Objectif
@objective(model, Max, sum(w[v]*x[v,1] for v in V))

println("==============================")
println("Formulation du modèle terminée")
println("==============================")
println("Nombre de variables : ", num_variables(model))
println("Objectif : ", objective_function(model))

# ----------------------------
# Optimisation en Branch-and-Cut
# ----------------------------

println("==============================")
println("Début de l'optimisation Branch-and-Cut")
println("==============================")

# 1. Configuration et enregistrement du callback
JuMP.set_optimizer_attribute(model, "OutputFlag", 1)
JuMP.set_optimizer_attribute(model, "LazyConstraints", 1) 

# --- CORRECTION DE LA SIGNATURE DU CALLBACK ---
# La fonction anonyme doit accepter deux arguments : cb_data et l'emplacement (cb_where)
# Le deuxième argument est ignoré dans notre appel à `component_separation_callback`
JuMP.set_optimizer_attribute(model, Gurobi.CallbackFunction(), (cb_data, cb_where) -> 
    component_separation_callback(cb_data, x, G, V, k; thr=0.5)
)
# ----------------------------------------------

# 2. Lancement de l'optimisation. Le solveur gère tout le processus B&C.
t1 = @elapsed optimize!(model)

# 3. Affichage des résultats
println("==============================")
println("Résultat final Branch-and-Cut")
println("==============================")
println("Statut : ", termination_status(model))
println("Valeur obj : ", objective_value(model))
println("Temps Branch-and-Cut : ", t1, " secondes")