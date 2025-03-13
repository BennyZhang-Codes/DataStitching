# push!(LOAD_PATH, "../src/")
using Documenter, HighOrderMRI

makedocs(
    sitename = "HighOrderMRI.jl",
    authors  = "Jinyuan Zhang, Zihao Zhang, Xiaoping Wu, and co-authors",
    # modules  = [HighOrderMRI],
    clean    = false,
    doctest  = false,
    pages = [
        "Introduction"  => "index.md",
        "API Reference" =>[
            "Reconstruction"  => "api/reconstruction.md",
            "Simulation"      => "api/simulation.md",
            "Synchronization" => "api/synchronization.md",
            "Plot function"   => "api/plt.md",
        ]
    ],
)


deploydocs(
    repo = "github.com/BennyZhang-Codes/HighOrderMRI.jl.git",
    devbranch = "main",  
    branch = "gh-pages",
    push_preview = true,
    versions = ["dev" => "main"]
)