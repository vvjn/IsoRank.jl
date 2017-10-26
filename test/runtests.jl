using IsoRank
using Base.Test

for k = 1:2:100
    G = sprand(k,k,0.4); G = G'+G; G = Int.(G - Diagonal(G) .> 0)
    R, res, L = isorank(G,G,maxiter=15,tol=1e-9,details=true,verbose=false)
    @test norm(L*vec(R) - res[1].*vec(R)) < 0.01
    b = rand(k,k);
    R, res, L = isorank(G,G,b,0.5,maxiter=25,tol=1e-9,details=true,verbose=false)
    @test norm(L*vec(R) - res[1].*vec(R)) < 0.01
end

function test1()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R, res, L = isorank(G,G,maxiter=15,tol=1e-9,details=true,verbose=false)
    norm(L*vec(R) .- res[1].*vec(R)) < 1e-4
end
@test test1()

function test2()
    G = sparse([1,2],[2,1],1,2,2)
    R,res,L = isorank(G, G, ones(2,2), 0.85, maxiter=20, tol=1e-5, details=true,verbose=false)
    vec(R) == [0.25,0.25,0.25,0.25] && norm(L*vec(R) - res[1].*vec(R)) == 0.0
end
@test test2()

function test3()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R,res,L = isorank(G, G, ones(size(G,1),size(G,1)), 0.85, maxiter=20, tol=1e-5, details=true,verbose=false)
    sum(vec(R) .== ones(9)./9)==0 && norm(L*vec(R) - res[1].*vec(R)) < 1e-10
end
@test test3()