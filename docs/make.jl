
using Documenter
using MarketModel

makedocs(sitename="MarketModel.jl",
    authors = "Richard Weinhold",
    pages = [
        "Introducion" => "index.md",
        ],
    );

deploydocs(
    repo = "github.com/richard-weinhold/MarketModel.git",
    devbranch = "construction",
    branch = "gh-pages"
)

