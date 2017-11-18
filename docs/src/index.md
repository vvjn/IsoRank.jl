# Introduction

[IsoRank.jl](https://github.com/vvjn/IsoRank.jl) is a Julia
implementation of IsoRank as described in "Global alignment of
multiple protein interaction networks with application to functional
orthology detection", Rohit Singh, Jinbo Xu, and Bonnie Berger
(2008). IsoRank.jl also contains a PageRank implementation. The greedy
network alignment method is also implemented here.

IsoRank calculates the topological similarity of all pairs of nodes
across two networks, with the assumption that a node is similar to
another node if the node's neighbors are similar to the other node's
neighbors. IsoRank can also use prior node similarity information
while calculating this topological node similarity measure.

The IsoRank matrix is calculated by creating the product graph of two
networks, and then performing PageRank on the product graph. PageRank
is done by using the power method to calculate the dominant
eigenvector of the modified adjacency matrix of the product
graph. Since IsoRank.jl doesn't explicitly build the product graph in
order to perform power iteration, it has much better time and space
complexity compared to other implementations of IsoRank. This
implementation of IsoRank runs in `O(K(|V|^2+|V||E|))`, where `|V|`
and`|E|` are the number of nodes and edges in the two networks, and
`K` is the total number of iterations required to converge under the
power method.

## Installation

IsoRank can be installed as follows. We also use the NetalignUtils to
read networks and so we install it as follows.

```julia
Pkg.clone("https://github.com/vvjn/IsoRank.jl")
```

## Example usage

We generate a scale-free network and create an IsoRank matrix between
the network and itself . We use a damping factor of 0.85 in order to
calculate a good IsoRank matrix using just network topology. We load
the `LightGraphs` package to generate networks.

```julia
using IsoRank, LightGraphs

g1 = erdos_renyi(200,0.1)
g2 = g1

G1 = adjacency_matrix(g1)
G2 = adjacency_matrix(g2)

R = isorank(G1, G2, 0.85)

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

The resulting alignment `f` describes a node mapping such that node
`i` in `g1` is mapped to node `j` in `g2` if `f[i] = j`. Thus, we can
create the aligned node pairs as follows.

``` julia
hcat(1:length(f), f[1:length(f)])
```

## Using node similarities

Assuming we have a matrix of prior node similarities, we can calculate
the IsoRank matrix while incorporating external information. We treat
the the node similarities as the personalization vector in PageRank.
Here, `b` is a matrix of node similarities. Here, we
equally weigh topological node similarity and prior node similarity by
setting the `alpha` variable to `0.5`. `alpha` must lie between `0.0`
and `1.0`. To give more weight to topological node similarity,
increase the `alpha` variable up to `1.0`.

```julia
b = rand(size(G1,1), size(G2,1))
b[sub2ind(size(b), 1:length(f), 1:length(f))] = 1.0

R = isorank(G1, G2, 0.5, b)
```

## Other parameters

Maximum number of iterations and error tolerance can be set as follows.

```julia
R = isorank(G1, G2, 0.5, b, maxiter=20, tol=1e-10)
```

We can extract the modified adjacency matrix, `L`, of the product
graph as follows. `vec(R)` is the dominant eigenvector and `res[1]` is
the corresponding eigenvalue of `L`.

```julia
R,res,L = isorank(G1, G2, 0.85, details=true)

println(norm(L * vec(R) - res[1] * vec(R),1))
```
