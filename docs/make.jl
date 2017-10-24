using Documenter, IsoRank

makedocs()

deploydocs(
           deps   = Deps.pip("mkdocs", "python-markdown-math"),           
           repo = "github.com/vvjn/IsoRank.jl.git",
           julia = "0.6"
)

# makedocs(
#     format = :html,
#     sitename = "IsoRank",
#     pages = [
#         "index.md"
#     ]
# )

