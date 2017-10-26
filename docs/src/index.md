# Introduction

[IsoRank.jl](https://github.com/vvjn/IsoRank.jl) is a Julia implementation of IsoRank as
described in "Global alignment of multiple protein interaction
networks with application to functional orthology detection", Rohit
Singh, Jinbo Xu, and Bonnie Berger (2008). IsoRank.jl also contains
a PageRank implementation.

The IsoRank matrix is calculated by creating the product graph of two
networks, and then performing PageRank on the product graph. PageRank
is done by using the power method to calculate the dominant
eigenvector of the modified adjacency matrix of the product
graph. Since IsoRank.jl doesn't explicitly build the product graph in
order to perform power iteration, it has much better time and space
complexity compared to other implementations of IsoRank.

The greedy network alignment method described in the paper is also implemented here.

## Installation

IsoRank can be installed as follows. We also use the NetalignUtils to
read networks and so we install it as follows.

```julia
Pkg.clone("https://github.com/vvjn/IsoRank.jl")
Pkg.clone("https://github.com/vvjn/NetalignMeasures.jl")
Pkg.clone("https://github.com/vvjn/NetalignUtils.jl")
```

## Example usage

We load an example network from the `examples/` directory and create
an IsoRank matrix between the network and itself. Unlike the original
paper which performs no damping when using network topology alone, we
give it a damping factor of 0.85 in order to calculate a good
IsoRank matrix using just network topology.

```julia
using NetalignUtils
using IsoRank

G1 = readgw("0Krogan_2007_high.gw").G
G2 = G1

R = isorank(G1, G2, damping=0.85)

R ./= maximum(R)
truemap = 1:size(G2,1)
randmap = randperm(size(G2,1))
println(sum(R[sub2ind(size(R),truemap,truemap)]))
println(sum(R[sub2ind(size(R),truemap,randmap)]))
```

Given the IsoRank matrix, we perform greedy alignment as follows.

``` julia
f = greedyalign(R)
```

Given the alignment `f`, we construct the aligned node pairs and
save the node pairs to file as follows.

``` julia
nodepairs = hcat(t1.nodes, t2.nodes[f])

writedlm("yeast_yeast.aln", nodepairs)
```

## Using node similarities

Assuming we have a matrix of node similarities, we can calculate the
IsoRank matrix while incorporating external information through node
similarities.  Here, `b` is a matrix of node similarities (but,
obviously, use meaningful node similarities instead of random values).

```julia
b = rand(size(G1,1), size(G2,1))

R = isorank(G1, G2, b, 0.5)
```

## Other parameters

Maximum number of iterations and error tolerance can be set as follows.

```julia
R = isorank(G1, G2, b, 0.5, maxiter=20, tol=1e-5)
```

We can extract the modified adjacency matrix, `L`, of the product graph as follows.
`vec(R)` is the dominant eigenvector and `res[1]` is the corresponding eigenvalue of `L`.

```julia
R,res,L = isorank(G1, G2, damping=0.85, details=true)

println(norm(L * vec(R) - res[1] * vec(R),1))
```