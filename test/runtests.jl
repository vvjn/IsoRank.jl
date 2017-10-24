using IsoRank
using Base.Test

function test1()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R, res, L = isorank(G,G,maxiter=15,tol=1e-9,details=true)
    norm(L*vec(R) .- res[1].*vec(R)) < 1e-4
end
@test test1()

function test2()
    G = sparse([1,2],[2,1],1,2,2)
    R,res,L = isorank(G, G, ones(2,2), 0.85, maxiter=20, tol=1e-5, details=true)
    vec(R) == [0.25,0.25,0.25,0.25] && norm(L*vec(R) - res[1].*vec(R)) == 0.0
end
@test test2()

function test3()
    G = sparse([1,2,3,1,3,2],[2,1,1,3,2,3],1,3,3)
    R,res,L = isorank(G, G, ones(size(G,1),size(G,1)), 0.85, maxiter=20, tol=1e-5, details=true)
    sum(vec(R) .== ones(9)./9)==0 && norm(L*vec(R) - res[1].*vec(R)) < 1e-10
end
@test test3()
