__precompile__()

module IsoRank

using IterativeSolvers
using LinearMaps

export isorank, kronlm

"""
Kronecker product of A and B  stored as a linearmap
so that you don't have to create the actual matrix.
This is much faster than creating the matrix like
the original paper does: O(|E|) instead of O(|E|^2)
for each step of the power iteration where |E| is the
average number of edges in the graphs
"""
function kronlm(T::Type,A,B)
    f = (y,x) -> begin
        y[:] = vec(B * reshape(x, size(B,2), size(A,2)) * A')
        y end
    fc = (y,x) -> begin
        y[:] = vec(B' * reshape(x, size(B,1), size(A,1)) * A)
        y end
    m = size(A,1) * size(B,1)
    n = size(A,2) * size(B,2)
    LinearMap{T}(f, fc, m, n)
end
kronlm(A,B) = kronlm(promote_type(eltype(A),eltype(B)),A,B)

"""
Creates the IsoRank matrix.
Rohit Singh, Jinbo Xu, and Bonnie Berger. (2008) Global alignment of
multiple protein interaction networks with application to functional
orthology detection, Proc. Natl. Acad. Sci. USA, 105:12763-12768.
Note: normalizes x with L_2 instead of L_1 like the paper does.
"""
function isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
                 b::AbstractMatrix, alpha::Real,
                 maxiter, tol)
    A = kronlm(Float64,G2,G1)
    d = 1.0 ./ (A * ones(Float64,size(A,2))) # rows of A sum to 1
    if alpha != 1.0 && b != nothing
        b = b ./ sum(abs,b) # b sums to 1
        L = LinearMap{Float64}((y,x) -> begin
                               At_mul_B!(y, A, d .* x)
                               y .= alpha .* y .+ (1-alpha) .* vec(b)
                               y
                               end, size(A,1), size(A,2))
    else
        L = LinearMap{Float64}((y,x) -> begin
                               At_mul_B!(y, A, d .* x)
                               y .= alpha .* y
                               y
                               end, size(A,1), size(A,2))
    end
    x = copy(vec(b)) #ones(Float64,size(L,2)) #rand(Float64,size(L,2))
    x ./= norm(x)
    res = powm!(L, x, log=true, verbose=true, tol=tol, maxiter=maxiter)
    reshape(x, size(G1,1), size(G2,1)), res, L
end

"""
Input is two adjacency matrices, a matrix of node similarities,
and weight between edge and node conservation
"""
function isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
                 b::AbstractMatrix, alpha::Real;
                 maxiter = 15, tol = eps(Float64) * size(G1,1) * size(G2,1))
    isorank(G1,G2,b,alpha,maxiter,tol)[1]
end

"""
If you don't have node similarities, you can still create
a decent IsoRank matrix by doing damping like PageRank does.
This is unlike the original paper that creates a very
bad IsoRank matrix when b = 0 or alpha = 1
"""
function isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC;
                 damping = 0.85,
                 maxiter = 15, tol = eps(Float64) * size(G1,1) * size(G2,1))
    b = ones(Float64,size(G1,1),size(G2,1))
    isorank(G1,G2,b,damping,maxiter,tol)[1]   
end

end
