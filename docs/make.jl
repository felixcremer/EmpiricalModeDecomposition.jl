using Documenter
using EmpiricalModeDecomposition

makedocs(
    sitename = "EmpiricalModeDecomposition",
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages=["Home" => "index.md",
           "Empirical Mode Decomposition" => "emd.md",
            "Library" => "api.md"],
    modules = [EmpiricalModeDecomposition]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<repository url>"
)
