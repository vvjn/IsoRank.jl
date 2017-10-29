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
