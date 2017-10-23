using IsoRank
using NetalignUtils
using Base.Test

G1 = readgw("0Krogan_2007_high.gw").G
G2 = G1
function test1(damping)
    b = ones(Float64,size(G1,1),size(G2,1))
    R, res, L = isorank(G1,G2,b,damping,maxiter=15,tol=1e-9,details=true)
    norm(L*vec(R) .- res[1].*vec(R)) < 1e-4
end
@test test1(0.4)
@test test1(0.8)

function test2()
    G = sparse([1,2],[2,1],1,2,2)
    R,res,L = isorank(G, G, ones(2,2), 0.85, maxiter=20, tol=1e-5, details=true)
    vec(R) == [0.25,0.25,0.25,0.25] && norm(L*vec(R) - res[1].*vec(R)) == 0.0
end
@test test2()

