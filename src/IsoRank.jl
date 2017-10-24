__precompile__()

module IsoRank

using LinearMaps

export isorank, kronlm, powermethod!

"""
    kronlm([::Type{T}], A, B)

Kronecker product of A and B  stored as a linear operator (from LinearMaps)
so that you don't have to create the actual matrix.
This is much faster than creating the matrix like
the original paper does: O(|E|) instead of O(|E|^2)
for each step of the power iteration where |E| is the
average number of edges in the graphs

# Arguments
- `A,B` : linear operators with multiply and transpose operations
- `T` : element type of the resulting linear operator  
"""
function kronlm(::Type{T},A,B) where {T}
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
    powermethod!(L, x; <keyword arguments>) -> radius, x, [log/history]

Performs power method in order to find the dominant eigenvector
of the linear operator L. Eigenvector is normalized w.r.t. L_1 norm.
Modifies initial eigenvector estimate x.

# Arguments
- `L` : linear operator
- `x` : initial estimate of the eigenvector, not necessarily normalized

# Keyword arguments
- `maxiter` : maximum # of iterations
- `tol=eps(Float64) * size(G1,1) * size(G2,1)` : error tolerance in L_1
- `log=true,verbose=true` : logging and printing
"""    
function powermethod!(A, x;
                      maxiter=15, tol=eps(Float64) * size(A,2),
                      log=true, verbose=true)
    T = eltype(A)
    x ./= norm(x,1)
    Ax = similar(x)
    verbose && println("Running power method, maxiter = $maxiter, tol = $tol")
    history = Tuple{Int,T,T,T}[]
    iter = 0
    radius = zero(T)
    lambda = zero(T)
    while iter <= maxiter
        A_mul_B!(Ax, A, x)
        lambda = norm(Ax,1)
        Ax ./= lambda # want |x|_1 = 1
        err = norm(Ax-x,1)
        radius = dot(x,Ax)/dot(x,x)
        verbose && @show iter,err,radius,lambda
        log && push!(history,(iter,err,radius,lambda))
        copy!(x, Ax)
        err < tol && break
        iter += 1
    end
    if log lambda,x,history else lambda,x end
end

"""
    isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
            b::AbstractMatrix, alpha::Real; <keyword arguments>)
    
Creates the IsoRank matrix. That is, finds the pagerank values
of the modified adjacency matrix of the product graph of G1 and G2.
b acts as node restarts allowing you to incorporate external information.    
(See Rohit Singh, Jinbo Xu, and Bonnie Berger. (2008) Global alignment of
multiple protein interaction networks with application to functional
orthology detection, Proc. Natl. Acad. Sci. USA, 105:12763-12768.)

# Arguments
- `G1,G2` : two adjacency matrices
- `b` : matrix of node similarities, not necessarily normalized
- `alpha`: weight between edge and node conservation

# Keyword arguments
- `details=false` : if true, returns (R,res,L) where R is the IsoRank
  matrix, res is the power method details structure, L is the linear
  operator that the power method finds the eigenvector of; if false,
  returns R
- See [`powermethod!`](@ref) for other keyword arguments   
"""
function isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
                 b::AbstractMatrix, alpha::Real;
                 details=false, args...)
    A = kronlm(Float64,G2,G1)
    d = 1.0 ./ (A * ones(Float64,size(A,2))) # rows of A sum to 1
    if alpha != 1.0
        bsum = norm(b,1)
        bsum==0.0 && error("b is a 0-vector")
        b = b ./ bsum # make b sum to 1
        L = LinearMap{Float64}((y,x) -> begin
                               At_mul_B!(y, A, d .* x)
                               y .= alpha .* y .+ (1.0-alpha) .* vec(b)
                               y
                               end, size(A,1), size(A,2))
        x = copy(vec(b))
    else
        L = LinearMap{Float64}((y,x) -> begin
                               At_mul_B!(y, A, d .* x)
                               y .= alpha .* y
                               y
                               end, size(A,1), size(A,2))
        x = ones(Float64,size(L,2)) #rand(Float64,size(L,2))
    end
    res = powermethod!(L, x; args...)
    R = reshape(x, size(G1,1), size(G2,1))
    if details
        R, res, L
    else
        R
    end
end

"""
    isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
            [damping=0.85]; <keyword arguments>)
    
If you don't have node similarities, you can still create
a decent IsoRank matrix by doing damping like PageRank does.
This is unlike the original paper that creates a
bad IsoRank matrix when b = 0 or alpha = 1.
"""
function isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC;
                 damping=0.85, args...)                 
    b = ones(Float64,size(G1,1),size(G2,1))
    isorank(G1,G2,b,damping; args...)   
end

end
