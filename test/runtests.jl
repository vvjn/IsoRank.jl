using IsoRank
using Base.Test

for k = 1:2:100
    G = sprand(k,k,0.4); G = G'+G; G = Int.(G - Diagonal(G) .> 0)
    R, res, L = isorank(G,G,maxiter=25,tol=1e-9,details=true,verbose=false)
    @test norm(L*vec(R) - res[1].*vec(R)) < 0.1
    b = rand(k,k);
    R, res, L = isorank(G,G,0.5,b,maxiter=25,tol=1e-9,details=true,verbose=false)
    @test norm(L*vec(R) - res[1].*vec(R)) < 0.1
end

function test1()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R, res, L = isorank(G,G,maxiter=25,tol=1e-9,details=true,verbose=false)
    norm(L*vec(R) .- res[1].*vec(R)) < 1e-2
end
@test test1()

function test2()
    G = sparse([1,2],[2,1],1,2,2)
    R,res,L = isorank(G, G, 0.85, ones(2,2), maxiter=25, tol=1e-9, details=true,verbose=false)
    vec(R) == [0.25,0.25,0.25,0.25] && norm(L*vec(R) - res[1].*vec(R)) == 0.0
end
@test test2()

function test3()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R,res,L = isorank(G, G, 0.85, ones(size(G,1),size(G,1)), maxiter=25, tol=1e-9, details=true,verbose=false)
    sum(vec(R) .== ones(9)./9)==0 && norm(L*vec(R) - res[1].*vec(R)) < 1e-2
end
@test test3()

function test4()
    diamond = sparse(vec([1 1 2 2 3]), vec([2 3 3 4 4]), 1.0,4,4)
    x = IsoRank.pagerank(diamond,verbose=false)
    isapprox(x, vec([0.128414 0.182991 0.260762 0.427833]), atol=0.01)
end
@test test4()

function test5()
    G1 = sprand(10,10,0.4); G1 = G1'+G1; G1 = Int.(G1 - Diagonal(G1) .> 0);
    G2 = sprand(15,15,0.4); G2 = G2'+G2; G2 = Int.(G2 - Diagonal(G2) .> 0);
    R, res, L = isorank(G1,G2,1.0,maxiter=25,tol=1e-9,details=true,verbose=false);
    A = Float64.(kron(G2,G1)); S = 1.0 ./ (A' * ones(Float64,size(A,2))); scale!(A,S);
    @test norm(A*vec(R) - res[1].*vec(R)) < 0.001
    els,evs = eigs(A);
    @test norm(A*evs[:,1] - L*evs[:,1]) < 0.001
end
test5()

    
