"""
Miscealenous functions for genetic genealogy.
Function to calculate an ethnicity estimate for MyHeritage matches files.
Functions to manipulate phylogenetic trees and to calculate TMRCAs.
"""
module GeneticGenealogy

export ethnicity, readtree, readtrees, print, writetree, calculateTMRCAs

include("ethnicity.jl")
include("tree.jl")

end # module
