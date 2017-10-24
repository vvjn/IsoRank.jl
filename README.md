# IsoRank

[![Build Status](https://travis-ci.org/vvjn/IsoRank.jl.svg?branch=master)](https://travis-ci.org/vvjn/IsoRank.jl) [![Coverage Status](https://coveralls.io/repos/vvjn/IsoRank.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/vvjn/IsoRank.jl?branch=master) [![codecov.io](http://codecov.io/github/vvjn/IsoRank.jl/coverage.svg?branch=master)](http://codecov.io/github/vvjn/IsoRank.jl?branch=master)

IsoRank.jl calculates the IsoRank matrix as described in "Global alignment of
multiple protein interaction networks with application to functional
orthology detection", Rohit Singh, Jinbo Xu, and Bonnie Berger (2008),
with much better space and time complexity. This is not code by the authors of
that paper.

The IsoRank matrix is calculated by creating the product graph of two
networks, and then performing PageRank on the product graph. PageRank
is calculated by performing power iteration to calculate the dominant
eigenvector of the modified adjacency matrix of the product
graph. Since IsoRank.jl doesn't explicitly build the product graph in
order to perform power iteration, it has much better time and space
complexity compared to other implementations of IsoRank.

## Installation

IsoRank can be installed as follows. We also use the NetalignUtils to
read networks.

```julia
Pkg.add("IsoRank")
Pkg.add("NetalignUtils")
```

## Example usage

We load an example network from the "test/" directory and create
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

truemap = 1:size(G2,1)
randmap = randperm(size(G2,1))
println(sum(R[sub2ind(size(R),truemap,truemap)]))
println(sum(R[sub2ind(size(R),truemap,randmap)]))
```

Assuming we have a matrix of node similarities, we can calculate
the IsoRank matrix using node similarities as follows, where `b` is
a matrix of node similarities.

```julia
b = rand(size(G1,1), size(G2,1))

R = isorank(G1, G2, b, 0.5)
```

Maximum number of iterations and error tolerance can be set as follows.

```julia
R = isorank(G1, G2, b, 0.5, maxiter=20, tol=1e-5)
```
