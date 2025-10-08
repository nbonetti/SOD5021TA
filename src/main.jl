using JuMP, Gurobi

# ----------------------
# Choix de l'instance et du modèle
# ----------------------
instance = "G_ex_papier.txt"
model = ""

# ----------------------
# Lire l'instance
# ----------------------
E_list, W = readWeightedGraph_paper(instance)
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
# Initialiser le modèle
# ----------------------


println("\n\n MODEL : ", model, "\n\n")


# ----------------------
# Résolution
# ----------------------
JuMP.optimize!(model)
println("Status = ", JuMP.termination_status(model))
println("Valeur optimale = ", JuMP.objective_value(model))

