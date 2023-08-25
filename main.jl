import Gmsh: gmsh
using SparseArrays
include("CreateMSH3D.jl")
include("FEM_3D_LIB.jl")

N_nodes, N_tetra, POS, TETRAHEDRA, TRIANGLES = CreateMSH3D("Cube")
k_dim = length(TETRAHEDRA[1,:])
weights, xs, ys, zs = get_tetra_quadrature_data(2)

ρ = 0
ε₀ = 8.854e-12
εᵣ = 1

BC = Dict(
    "D"=>Dict("flag"=>"Neumann", "val"=>0), 
    "U"=>Dict("flag"=>"Neumann", "val"=>0), 
    "F"=>Dict("flag"=>"Neumann", "val"=>0),  
    "R"=>Dict("flag"=>"Dirichlet", "val"=>10),  
    "B"=>Dict("flag"=>"Neumann", "val"=>0),  
    "L"=>Dict("flag"=>"Dirichlet", "val"=>0)
)

Dirichlet_nodes_vector = UInt64[]
Dirichlet_values_vector = UInt64[]
for dir ∈ ["D", "U", "F", "R", "B", "L"]
    if BC[dir]["flag"] == "Dirichlet"
        for n ∈ TRIANGLES[dir]
            if n ∉ Dirichlet_nodes_vector
                append!(Dirichlet_nodes_vector, n)
                append!(Dirichlet_values_vector, BC[dir]["val"])
            end   
        end
    end
end

I = zeros(N_tetra*(k_dim)^2)
J = zeros(N_tetra*(k_dim)^2)
V = zeros(N_tetra*(k_dim)^2)
F = zeros(N_nodes)
for el ∈ 1:N_tetra
    nodes_numbers = TETRAHEDRA[el,:]
    X = POS[nodes_numbers,:]'
    g(ξ₁, ξ₂, ξ₃) = ∇N∇Nᵀ_func(ξ₁, ξ₂, ξ₃, X)
    h(ξ₁, ξ₂, ξ₃) = N_func(ξ₁, ξ₂, ξ₃, X)
    k = tetra_quadrature(g, (k_dim,k_dim), weights, xs, ys, zs)
    f = tetra_quadrature(h, k_dim, weights, xs, ys, zs) * ρ / (ε₀ * εᵣ)
    for i = 1:4
        I[16*el-(15-4*(i-1)) : 16*el-(12-4*(i-1))] .+= nodes_numbers
        J[16*el-(15-4*(i-1)) : 16*el-(12-4*(i-1))] .+= nodes_numbers[i]
    end
    V[16*(el-1)+1:16*el] .+= reshape(k,(16,))
    F[nodes_numbers] .+= f
end
K = sparse(I, J, V, N_nodes, N_nodes, +)

for dir ∈ ["D", "U", "F", "R", "B", "L"]
    if BC[dir]["flag"] == "Neumann"
        for i_tri ∈ 1:size(TRIANGLES[dir],1)
            nodes_numbers = TRIANGLES[dir][i_tri,:]
            M = POS[nodes_numbers,:]'
            F[nodes_numbers] .-= 2 * BC[dir]["val"] * TriangleArea3D(M) * [1; 1; 1] / 6 
        end
    end
end

for (i, dir_node_num) ∈ enumerate(Dirichlet_nodes_vector)
    K[dir_node_num,:] .= 0
    K[dir_node_num,dir_node_num] = 1
    F[dir_node_num] = Dirichlet_values_vector[i]
end

dropzeros!(K)

ϕ = K \ F

err = zeros(N_nodes)
for i ∈ 1:N_nodes
    err[i] = ϕ[i] - (5*POS[i,1]+5)
end

println(findmax(err))
