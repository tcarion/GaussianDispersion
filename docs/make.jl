using GaussianDispersion
using Documenter

DocMeta.setdocmeta!(GaussianDispersion, :DocTestSetup, :(using GaussianDispersion); recursive=true)

makedocs(;
    modules=[GaussianDispersion],
    authors="tcarion <tristan.carion@gmail.com> and contributors",
    repo="https://github.com/tcarion/GaussianDispersion.jl/blob/{commit}{path}#{line}",
    sitename="GaussianDispersion.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tcarion.github.io/GaussianDispersion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tcarion/GaussianDispersion.jl",
    devbranch="main",
)
