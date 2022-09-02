using Equil
using Documenter

DocMeta.setdocmeta!(Equil, :DocTestSetup, :(using Equil); recursive=true)

makedocs(;
    modules=[Equil],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/Equil.jl/blob/{commit}{path}#{line}",
    sitename="Equil.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/Equil.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/Equil.jl",
    devbranch="main",
)
