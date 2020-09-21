cd(@__DIR__)
using Pkg
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
Pkg.activate(@__DIR__)
CI && Pkg.instantiate()
using ClimateBase
using Documenter, DimensionalData, NCDatasets, StatsBase

# %%
# using DocumenterTools: Themes
# download the themes
# for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
#     download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
# end
# # create the themes
# for w in ("light", "dark")
#     header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
#     theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
#     write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
# end
# # compile the themes
# Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
# Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))


# %% Build the actual docs
makedocs(modules = [ClimateBase, DimensionalData],
    sitename= "ClimateBase",
    authors = "George Datseris, Phillippe Roy, and contributors.",
    doctest = false,
    format = Documenter.HTML(
        prettyurls = CI,
        # assets = [
        #     "assets/logo.ico",
        #     asset("https://fonts.googleapis.com/css?family=Quicksand|Montserrat|Source+Code+Pro|Lora&display=swap", class=:css),
        #     ],
        ),
    pages = [
        "Introduction" => "index.md",
        ],
)

if CI
    deploydocs(
        repo = "github.com/JuliaClimate/ClimateBase.jl.git",
        target = "build",
        push_preview = true
    )
end
