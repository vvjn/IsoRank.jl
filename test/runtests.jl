using IsoRank
using NetalignUtils
using Base.Test

G1 = readnetwork("0Krogan_2007_high.gw").G
G2 = readnetwork("0Krogan_2007_high.gw").G

for damping in [0.4, 0.8]
    b = ones(Float64,size(G1,1),size(G2,1))
    R, res, L = isorank(G1,G2,b,damping,15,1e-9)
    @test abs(res[1] - damping) < 0.01
end


