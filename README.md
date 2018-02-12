# IsoRank

[![Build Status](https://travis-ci.org/vvjn/IsoRank.jl.svg?branch=master)](https://travis-ci.org/vvjn/IsoRank.jl) [![codecov.io](http://codecov.io/github/vvjn/IsoRank.jl/coverage.svg?branch=master)](http://codecov.io/github/vvjn/IsoRank.jl?branch=master) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://vvjn.github.io/IsoRank.jl/latest)

IsoRank.jl is a Julia implementation of IsoRank as
described in "Global alignment of multiple protein interaction
networks with application to functional orthology detection", Rohit
Singh, Jinbo Xu, and Bonnie Berger (2008). IsoRank.jl also contains
a PageRank implementation. The greedy network alignment method
is also implemented here.

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

IsoRank can be installed as follows.

```julia
Pkg.add("IsoRank")
```

# Documentation

Available [here](https://vvjn.github.io/IsoRank.jl/latest).
