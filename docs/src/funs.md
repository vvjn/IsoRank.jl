# Functions


```@meta
CurrentModule = IsoRank
```

```@docs
isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,
         b::AbstractMatrix, alpha::Real;
         details=false, args...)
```

```@docs
isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC;
         damping=0.85, args...)                 
```

```@docs
powermethod!(A, x; maxiter=15, tol=eps(Float64) * size(A,2),
                    log=true, verbose=true)
```

```@docs
kronlm(::Type{T},A,B) where {T}
```
