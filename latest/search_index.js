var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "IsoRank.jl is a Julia implementation of IsoRank as described in \"Global alignment of multiple protein interaction networks with application to functional orthology detection\", Rohit Singh, Jinbo Xu, and Bonnie Berger (2008). IsoRank.jl also contains a PageRank implementation.The IsoRank matrix is calculated by creating the product graph of two networks, and then performing PageRank on the product graph. PageRank is done by using the power method to calculate the dominant eigenvector of the modified adjacency matrix of the product graph. Since IsoRank.jl doesn't explicitly build the product graph in order to perform power iteration, it has much better time and space complexity compared to other implementations of IsoRank.The greedy network alignment method described in the paper is also implemented here."
},

{
    "location": "index.html#Installation-1",
    "page": "Introduction",
    "title": "Installation",
    "category": "section",
    "text": "IsoRank can be installed as follows. We also use the NetalignUtils to read networks and so we install it as follows.Pkg.clone(\"https://github.com/vvjn/IsoRank.jl\")\nPkg.clone(\"https://github.com/vvjn/NetalignMeasures.jl\")\nPkg.clone(\"https://github.com/vvjn/NetalignUtils.jl\")"
},

{
    "location": "index.html#Example-usage-1",
    "page": "Introduction",
    "title": "Example usage",
    "category": "section",
    "text": "We load an example network from the examples/ directory and create an IsoRank matrix between the network and itself. Unlike the original paper which performs no damping when using network topology alone, we give it a damping factor of 0.85 in order to calculate a good IsoRank matrix using just network topology.using NetalignUtils\nusing IsoRank\n\nG1 = readgw(\"0Krogan_2007_high.gw\").G\nG2 = G1\n\nR = isorank(G1, G2, damping=0.85)\n\nR ./= maximum(R)\ntruemap = 1:size(G2,1)\nrandmap = randperm(size(G2,1))\nprintln(sum(R[sub2ind(size(R),truemap,truemap)]))\nprintln(sum(R[sub2ind(size(R),truemap,randmap)]))Given the IsoRank matrix, we perform greedy alignment as follows.f = greedyalign(R)Given the alignment f, we construct the aligned node pairs and save the node pairs to file as follows.nodepairs = hcat(t1.nodes, t2.nodes[f])\n\nwritedlm(\"yeast_yeast.aln\", nodepairs)"
},

{
    "location": "index.html#Using-node-similarities-1",
    "page": "Introduction",
    "title": "Using node similarities",
    "category": "section",
    "text": "Assuming we have a matrix of node similarities, we can calculate the IsoRank matrix while incorporating external information through node similarities.  Here, b is a matrix of node similarities (but, obviously, use meaningful node similarities instead of random values).b = rand(size(G1,1), size(G2,1))\n\nR = isorank(G1, G2, b, 0.5)"
},

{
    "location": "index.html#Other-parameters-1",
    "page": "Introduction",
    "title": "Other parameters",
    "category": "section",
    "text": "Maximum number of iterations and error tolerance can be set as follows.R = isorank(G1, G2, b, 0.5, maxiter=20, tol=1e-5)We can extract the modified adjacency matrix, L, of the product graph as follows. vec(R) is the dominant eigenvector and res[1] is the corresponding eigenvalue of L.R,res,L = isorank(G1, G2, damping=0.85, details=true)\n\nprintln(norm(L * vec(R) - res[1] * vec(R),1))"
},

{
    "location": "funs.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "funs.html#IsoRank.isorank-Tuple{SparseMatrixCSC,SparseMatrixCSC,AbstractArray{T,2} where T,Real}",
    "page": "Functions",
    "title": "IsoRank.isorank",
    "category": "Method",
    "text": "isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,\n        b::AbstractMatrix, alpha::Real; <keyword arguments>)\n\nCreates the IsoRank matrix. That is, finds the PageRank values of the modified adjacency matrix of the product graph of G1 and G2. b acts as node restarts allowing you to incorporate external information.     (See Rohit Singh, Jinbo Xu, and Bonnie Berger. (2008) Global alignment of multiple protein interaction networks with application to functional orthology detection, Proc. Natl. Acad. Sci. USA, 105:12763-12768.)\n\nArguments\n\nG1,G2 : two adjacency matrices\nb : matrix of node similarities between G1 and G2, not necessarily normalized\nalpha: weight between edge and node conservation\n\nKeyword arguments\n\ndetails=false : If true, returns (R,res,L) where R is the IsoRank matrix, res is the power method detailed results structure, L is the linear operator that the power method finds the eigenvector of; if false, returns R\nSee powermethod! for other keyword arguments   \n\n\n\n"
},

{
    "location": "funs.html#IsoRank.isorank-Tuple{SparseMatrixCSC,SparseMatrixCSC}",
    "page": "Functions",
    "title": "IsoRank.isorank",
    "category": "Method",
    "text": "isorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,\n        [damping=0.85]; <keyword arguments>)\n\nIf you don't have node similarities, you can still create a decent IsoRank matrix by doing damping like PageRank does. This is unlike the original paper that creates a bad IsoRank matrix when b = 0 or alpha = 1.\n\n\n\n"
},

{
    "location": "funs.html#IsoRank.greedyalign",
    "page": "Functions",
    "title": "IsoRank.greedyalign",
    "category": "Function",
    "text": "greedyalign(R::AbstractMatrix,seeds=Vector{Tuple{Int,Int}}();\n            maxiter=size(R,1)-length(seeds))\n\nGiven matrix R, find a alignment using the greedy method.\n\nArguments\n\nR : m x n matrix\nseeds : Initial seed node pairs\nmaxiter : Stop after aligning maxiter node pairs\n\n\n\n"
},

{
    "location": "funs.html#IsoRank.pagerank",
    "page": "Functions",
    "title": "IsoRank.pagerank",
    "category": "Function",
    "text": " pagerank(A, damping=0.85, p = fill(1/size(A,2),size(A,2));\n          <keyword args>) -> x [, res, L]\n\nCreates PageRank vector.\n\nArguments\n\nA : Adjacency matrix of the graph. A[u,v] = true if u -> v\ndamping : Damping\np : Initial probability vector\n\nKeyword arguments\n\ndetails=false : If true, returns (x,res,L) where x is the PageRank vector, res is the power method detailed results structure, L is the linear operator that the power method finds the eigenvector of; if false, returns x.\nSee powermethod! for other keyword arguments   \n\n\n\n"
},

{
    "location": "funs.html#IsoRank.powermethod!-Tuple{Any,Any}",
    "page": "Functions",
    "title": "IsoRank.powermethod!",
    "category": "Method",
    "text": "powermethod!(A, x; <keyword arguments>) -> radius, x, [log/history]\n\nPerforms power method in order to find the dominant eigenvector of the linear operator A. Eigenvector is normalized w.r.t. L_1 norm. Modifies initial eigenvector estimate x.\n\nArguments\n\nA : linear operator\nx : initial estimate of the eigenvector, not necessarily normalized\n\nKeyword arguments\n\nmaxiter : maximum # of iterations\ntol=eps(Float64) * size(A,2) : error tolerance in L_1 norm\nlog=true,verbose=true : logging and printing\n\n\n\n"
},

{
    "location": "funs.html#IsoRank.kronlm-Union{Tuple{Type{T},Any,Any}, Tuple{T}} where T",
    "page": "Functions",
    "title": "IsoRank.kronlm",
    "category": "Method",
    "text": "kronlm([T], A, B)\n\nKronecker product of A and B, stored as a linear operator (from LinearMaps.jl) so that you don't have to create the actual matrix. This is much faster than creating the matrix like the original paper does: O(|E|) instead of O(|E|^2) for each step of the power iteration where |E| is the average number of edges in the graphs\n\nArguments\n\nA,B : linear operators with multiply and transpose operations\nT : element type of the resulting linear operator \n\n\n\n"
},

{
    "location": "funs.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = IsoRankisorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC,\n         b::AbstractMatrix, alpha::Real;\n         details=false, args...)\nisorank(G1::SparseMatrixCSC, G2::SparseMatrixCSC;\n         damping=0.85, args...)\ngreedyalign\npagerank\npowermethod!(A, x; maxiter=15, tol=eps(Float64) * size(A,2),\n                    log=true, verbose=true)\nkronlm(::Type{T},A,B) where {T}"
},

]}
