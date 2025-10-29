using JuMP, Gurobi
include("flow.jl")
include("asymmetric_flow.jl")
include("cut_based_connectivity.jl")
include("cut_based_cross.jl")
include("utils.jl")

# ----------------------
# Choix de l'instance et du modèle
# ----------------------
instance = "G_ex_papier.txt"
model = "flow"

# ----------------------
# Lire l'instance
# ----------------------
instance_path = "instances/" * instance
E_list, W = readWeightedGraph_paper(instance_path)
println("===== Recuperation de l'instance =====")
println("Instance lue ", instance)
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

if model == "flow"
    run_flow_model(E_list, W, Wtot, n, V, k)
elseif model == "asymmetric_flow"
    run_asymmetric_flow_model(E_list, W, Wtot, n, V, k)
elseif model == "cut_based_connectivity"
    run_cut_based_connectivity_model(E_list, W, n, V, k)
elseif model == "cut_based_cross"
    run_cut_based_cross_model(E_list, W, n, V, k)
end
