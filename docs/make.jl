using Documenter
using NLCE

makedocs(
    sitename = "NLCE",
    format = Documenter.HTML(),
    modules = [NLCE]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
