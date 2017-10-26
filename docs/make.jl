using Documenter, IsoRank

makedocs(
           format = :html,
           sitename = "Isorank",
           pages = [
                    "index.md", "Functions" => "funs.md"
           ]
       )

deploydocs(
           repo = "github.com/vvjn/IsoRank.jl.git",
           target = "build",
           deps   = nothing,
           make   = nothing,
           julia = "0.6"
)
