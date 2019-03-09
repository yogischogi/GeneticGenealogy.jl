import Base: print, replace!, show

# Indentation for output.
const INDENT_STEP = 4

abstract type TreeElement end

mutable struct Sample <: TreeElement
    fields::Array{String}
end

mutable struct Clade <: TreeElement
    name::String
    fields::Array{String}
    samples::Array{Sample}
    subclades::Array{Clade}
end

function Clade(fields::Array{<:AbstractString})
    name = strip(fields[1])
    Clade(name, fields, [], [])
end


function readtree(filename::AbstractString)
    # Read tree from file in text format.
    lines = String[]
    open(filename) do f
        for line in eachline(f)
            if !startswith(line, "//") && !all(isspace, line)
                push!(lines, line)
            end
        end
    end
    # Build tree.
    if !isempty(lines)
        root = parsetree!(lines)
    else
        error("Empty file.")
    end
end


function readtrees(filenames::Array{<:AbstractString})
    if isempty(filenames)
        error("No filenames provided.")
    end
    root = readtree(filenames[1])
    if length(filenames) > 1
        for filename in filenames[2:end]
            tree = readtree(filename)
            root = mergetrees!(root, tree)
        end
    end
    return root
end


function parsetree!(lines::Array{String})
    # Create top node.
    line = popfirst!(lines)
    indent = indentation(line)
    fields = split(strip(line), ",")
    node = Clade(fields)

    # Parse samples and subclades
    while !isempty(lines)
        if indentation(lines[1]) <= indent
            break
        end
        content = strip(lines[1])
        if startswith(content, "id:")
            fields = split(content, ",")
            push!(node.samples, Sample(fields))
            popfirst!(lines)
        else
            subclade = parsetree!(lines)
            push!(node.subclades, subclade)
        end
    end    
    return node
end

function Base.show(io::IO, node::Clade)
    print(io, node, 0)
end

function print(io::IO, node::Clade, indent::Int)
    subindent = indent + INDENT_STEP
    printspaces(io, indent)
    printfields(io, node.fields)
    println(io)
    for sample in node.samples
        printspaces(io, subindent)
        printfields(io, sample.fields)
        println(io)
    end
    for subclade in node.subclades
        print(io, subclade, subindent)
    end
end

function writetree(filename::String, tree::Clade)
    buffer = IOBuffer()
    print(buffer, tree, 0)
    open(filename, "w") do f
        write(f, String(take!(buffer)))
    end
end


function mergetrees!(tree::Clade, subtree::Clade)
    for i in eachindex(tree.subclades)
        if tree.subclades[i].name == subtree.name
            tree.subclades[i] = subtree
            break
        else
            mergetrees!(tree.subclades[i], subtree)
        end
    end
    return tree
end


function replace!(tree::Clade, pat_f::Pair)
    for i in eachindex(tree.fields)
        tree.fields[i] = replace(tree.fields[i], pat_f)
    end
    for sample in tree.samples
        for i in eachindex(sample.fields)
            sample.fields[i] = replace(sample.fields[i], pat_f)
        end
    end
    for subclade in tree.subclades
        replace!(subclade, pat_f)
    end
end

function deletefields!(tree::Clade, keys::Array{<:AbstractString})
    tree.fields = deletefields(tree.fields, keys)
    for sample in tree.samples
        sample.fields = deletefields(sample.fields, keys)
    end
    for subclade in tree.subclades
        deletefields!(subclade, keys)
    end
end


function deletefields(fields:: Array{<:AbstractString}, keys::Array{<:AbstractString})
    newfields = []
    for field in fields
        isclean = true
        for key in keys
            if occursin(key, field)
                isclean = false
                break
            end
        end
        if isclean
            push!(newfields, field)
        end
    end
    return newfields
end


# Needs Phyloage to be installed.
# https://github.com/yogischogi/phyloage
function calculateTMRCAs(tree::Clade; cal::Float64=140.0, offset::Float64=60.0)
    caltree = deepcopy(tree)
    replace!(caltree, "snps:" => "STR-Count:")
    
    # Call Phyloage with temporary files to calculate results.
    infile, inio = mktemp(pwd())
    outfile, outio = mktemp(pwd())
    close(inio)
    close(outio)
    writetree(infile, caltree)
    run(`phyloage -treein=$infile -treeout=$outfile -offset=$offset -cal=$cal`)
    resulttree = readtree(outfile)
    rm(infile)
    rm(outfile)

    # Remove clutter from result tree.
    deletefields!(resulttree, ["STR-Count:", "STRs Downstream:", "formed:"])
    return resulttree
end

# Count white spaces.
function indentation(line::String)
    indent = 0
    for char in line
        isspace(char) ? indent += 1 : break
    end
    return indent
end

function printspaces(io::IO, count::Int)
    for i = 1:count
        print(io, " ")
    end
end

function printfields(io::IO, fields::Array{<:AbstractString})
    if isempty(fields)
        return
    end
    for i in 1:length(fields)-1
        print(io, fields[i])
        print(io, ",")
    end
    print(io, fields[end])
end






