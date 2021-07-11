# Functions to calculate an ethnicity estimate.

using DataFrames
using DelimitedFiles
using Printf

const DEBUG = true

const CHROMOSOMES = 22

# Length of the chromosome segments that are used for mapping.
const SEGMENT_LENGTH = 1000000

# Matching segments that are larger than MAX_MATCH_LENGTH are
# excluded from the evaluation because they indicate close family
# members. We do not want them to obfuscate the results..
const MAX_MATCH_LENGTH = 20000000

# Sizes of the different chromosomes. We use larger numbers here
# because chromosome sizes can change with different genome reference versions.
const CHROMOSOME_SIZES = [250000000, 250000000, 200000000, 200000000, 190000000,
                   180000000, 160000000, 150000000, 150000000, 140000000,
                   140000000, 140000000, 120000000, 110000000, 110000000,
                   100000000,  90000000,  90000000,  60000000,  70000000,
                    50000000,  60000000]

"""
    ethnicity(matches, segments, [parent_matches,] birth_country="", excludes=["USA", "Canada", "Australia"])

Calculate an ethnicity estimate from a MyHeritage "DNA matches
list" and the corresponding "shared DNA segments" list.
If a birth country is provided the calculation will be more accurate.
You can specify a number of countries that are excluded from the calculation
because of massive recent migrations.

Important: The CSV files must be in MyHeritage's English file format!

# Examples
```
julia> ethnicity("matches.csv", "shared_segments.csv")
Total: 1311 segments
        Germany  465  35.5%
        Denmark  336  25.6%
    Netherlands  206  15.7%

julia> ethnicity("matches.csv", "shared_segments.csv", birth_country="Germany")
Total: 1124 segments
        Germany  710  63.2%
        Denmark  164  14.6%
    Netherlands  106   9.4%
```
"""
function ethnicity end


function ethnicity(matches_file::AbstractString,
                   segments_file::AbstractString;
                   birth_country = "",
                   excludes=["USA", "Canada", "Australia"])
    local matches
    local segments

    # Read files.
    try
        matches = open(matches_file) do f
            readdlm(f, ',', String)
        end
        segments = open(segments_file) do f
            readdlm(f, ',', String)
        end
    catch ex
        error("Could not read file. ", ex)
    end

    # Convert matches and segments to DataFrames.
    try
        matches = matches_to_dataframe(matches)
        segments = segments_to_dataframe(segments)
    catch ex
        error("Could not interpret files. Wrong format? ", ex)
    end

    # Map matches' segments onto DNA.
    matches_segments = innerjoin(matches, segments, on = :Name)
    dna = DNA()
    map_countries!(dna, matches_segments, birth_country, excludes)

    # Evaluate ethnicities.
    eths, segments_total = ethnicities(dna)
    print_ethnicities(eths, segments_total)
end


function ethnicity(matches_file::AbstractString,
                   segments_file::AbstractString,
                   parent_matches_file::AbstractString;
                   birth_country = "",
                   excludes=["USA", "Canada", "Australia"])
    local matches
    local segments
    local parent

    # Read files.
    try
        matches = open(matches_file) do f
            readdlm(f, ',', String)
        end
        segments = open(segments_file) do f
            readdlm(f, ',', String)
        end
        parent = open(parent_matches_file) do f
            readdlm(f, ',', String)
        end
    catch ex
        error("Could not read file. ", ex)
    end

    # Convert matches and segments to DataFrames.
    try
        matches = matches_to_dataframe(matches)
        segments = segments_to_dataframe(segments)
        parent = matches_to_dataframe(parent)
    catch ex
        error("Could not interpret files. Wrong format? ", ex)
    end

    # Phase matches.
    # Matches that match child and parent.
    child_parent = innerjoin(matches, parent, on = :Name, makeunique=true)
    # Matches that only match the second parent.
    child_parent2 = antijoin(matches, child_parent, on = :Name, makeunique=true)
    # Remove matches that match both parents from the child-parent matches.
    both, child_parent = remove_false_matches(child_parent, parent)

    # Add segments to the matches lists.
    child_parent_segments = innerjoin(child_parent, segments, on = :Name)
    child_parent2_segments = innerjoin(child_parent2, segments, on = :Name)
    both_segments = innerjoin(both, segments, on = :Name)

    # Map matches' segments onto DNA.
    parent_dna = DNA()
    parent2_dna = DNA()
    both_dna = DNA()
    map_countries!(parent_dna, child_parent_segments, birth_country, excludes)
    map_countries!(parent2_dna, child_parent2_segments, birth_country, excludes)
    map_countries!(both_dna, both_segments, birth_country, excludes)

    # Evaluate ethnicities.
    parent_eths, segments_total_1 = ethnicities(parent_dna)
    parent2_eths, segments_total_2 = ethnicities(parent2_dna)
    both_eths, segments_total_3 = ethnicities(both_dna)

    println("\nEthnicities you inherited from your parent.")
    print_ethnicities(parent_eths, segments_total_1)

    println("\nEthnicities you inherited from your other parent.")
    print_ethnicities(parent2_eths, segments_total_2)

    println("\nEthnicities that contributed to you and both of your parents.")
    print_ethnicities(both_eths, segments_total_3)
end

"""
DNA is modeled as an array of chromosomes.
Each chromosome consists of a number of segments represented
by a dictionary that maps country names to the number of matches
from each country.
"""
mutable struct DNA
    chromosomes::Array{Array{Dict{String, Integer}}}
end

function DNA()
    chromosomes = Array{Array{Dict{String, Integer}}}(undef, CHROMOSOMES)
    for i in 1:CHROMOSOMES
        nSegs = div(CHROMOSOME_SIZES[i] + SEGMENT_LENGTH - 1, SEGMENT_LENGTH)
        segments = Array{Dict{String, Integer}}(undef, nSegs)
        for j in 1:nSegs
            segments[j] = Dict{String, Integer}()
        end
        chromosomes[i] = segments
    end
    DNA(chromosomes)
end

"""
    map_countries!(dna::DNA, df::DataFrame)

Map country information from matches onto the DNA chromosomes.
The data frame "segments" must contain the following columns:
:Country, :Chromosome, :StartLocation, :EndLocation.
excludes contains a list of countries that should be excluded from mapping.
"""
function map_countries!(dna::DNA, segments::DataFrame, birth_country::String, excludes::Array{String})
    # Map matches' segments onto the DNA.
    nrows, _ = size(segments)
    for i in 1: nrows
        # Ignore close relatives.
        if segments[i, :EndLocation] - segments[i, :StartLocation] > MAX_MATCH_LENGTH
            continue
        end

        # Map matching segment onto the DNA.
        a = segments[i, :StartLocation] + SEGMENT_LENGTH
        b = segments[i, :EndLocation]
        c = segments[i, :Chromosome]
        country = segments[i, :Country]
        # Ignore countries that should be excluded.
        if in(country, excludes)
            continue
        end
        for j in a:SEGMENT_LENGTH:b
            seg = div(j, SEGMENT_LENGTH)
            if haskey(dna.chromosomes[c][seg], country)
                dna.chromosomes[c][seg][country] += 1
            else
                dna.chromosomes[c][seg][country] = 1
            end
        end
    end
    # Map birth country onto DNA segments that are already populated.
    if birth_country == "" || in(birth_country, excludes)
        return
    end
    for i in eachindex(dna.chromosomes)
        for j in eachindex(dna.chromosomes[i])
            # Check if DNA segment is populated.
            if length(keys(dna.chromosomes[i][j])) == 0
                continue
            end
            # Add birth country.
            if haskey(dna.chromosomes[i][j], birth_country)
                dna.chromosomes[i][j][birth_country] += 1
            else
                dna.chromosomes[i][j][birth_country] = 1
            end
        end
    end
end

"""
    ethnicities(dna::DNA)

Calculate ethnicities for a given DNA sample.

Returns an Array{Array{String}} and the total number of
segments that were evaluated. For each chromosome an array
exists that consists of segment entries. A segment entry is
a String that holds the most likely country origin of the
segment.
"""
function ethnicities(dna::DNA)
    # Create an empty country entry for each DNA segment.
    ethSegments = Array{Array{String}}(undef, CHROMOSOMES)
    for i in 1:CHROMOSOMES
        nSegs = length(dna.chromosomes[i])
        chromosome = Array{String}(undef, nSegs)
        for j in 1:nSegs
            chromosome[j] = ""
        end
        ethSegments[i] = chromosome
    end

    # Determine ethnicity for each segment.
    for i in 1:CHROMOSOMES
        for j in eachindex(dna.chromosomes[i])
            ethSegments[i][j] = most_common_ethnicity(dna.chromosomes[i][j])
        end
    end

    # Count segment ethnicities.
    ethnicities = Dict{String, Integer}()
    # Total number of evaluated segments.
    total = 0
    for i in 1:CHROMOSOMES
        for j in eachindex(ethSegments[i])
            if ethSegments[i][j] != ""
                total += 1
                if haskey(ethnicities, ethSegments[i][j])
                    ethnicities[ethSegments[i][j]] += 1
                else
                    ethnicities[ethSegments[i][j]] = 1
                end
            end
        end
    end
    return ethnicities, total
end

"""
    print_ethnicities(eth::Dict{String, Integer}, segments_total::Integer)

Print the result of the ethnicities function.
"""
function print_ethnicities(eth::Dict{String, Integer}, segments_total::Integer)
    # Sort results.
    results = []
    for (country, value) in eth
        percentage = value / segments_total * 100
        push!(results, (country, value, percentage))
    end
    sort!(results, by = x -> x[2], rev = true)

    # Output results.
    println("Total: ", segments_total, " segments")
    for i in results
        @printf("%15s %4d %5.1f%%\n", i[1], i[2], i[3])
    end
end

"""
    matches_to_dataframe(matches::Array{String, 2})

Converts matches to a valid DataFrame. MyHeritage CSV files
sometimes contain an invalid CSV format and can not be read by
the standard Julia CSV package. So we must do this manually.
"""
function matches_to_dataframe(matches::Array{String, 2})
    df = DataFrame(Name = String[], Country = String[], Total_cM_shared = Float64[])

    # Copy matches rows into DataFrame.
    names = Dict{String, Int64}()
    rows, _ = size(matches)
    for i = 2:rows
        country = matches[i, 4]
        name = matches[i, 2]
        matches_cM = matches[i, 10]
        if country == "" || name == "" || matches_cM == ""
            continue
        end

        # Ensure that names are unique.
        if !haskey(names, name)
            names[name] = 1
        else
            names[name] += 1
            name = name * string(names[name])
        end

        # Convert cM to Float.
        cM = parse(Float64, replace(matches_cM, "," => ""))
        push!(df, [name, country, cM])
    end
    return df
end

"""
    segments_to_dataframe(segments::Array{String, 2})

Converts matches to a valid DataFrame. MyHeritage CSV files
sometimes contain an invalid CSV format and no valid key column.
They can not always be read by the standard Julia CSV package.
So we must do this manually.
"""
function segments_to_dataframe(segments::Array{String, 2})
    df = DataFrame(Name = String[], Chromosome = Integer[], StartLocation = Integer[], EndLocation = Integer[])

    # Copy segments rows into DataFrame.
    names = Dict{String, Int64}()
    prev_name = ""
    block_name = ""
    rows, _ = size(segments)
    for i = 2:rows
        name = segments[i, 3]
        chromosome = segments[i, 4]
        startloc = segments[i, 5]
        endloc = segments[i, 6]
        if name == "" || chromosome == "" || startloc == "" || endloc == ""
            continue
        end
        
        # Make sure that different persons get different names.
        if !haskey(names, name)
            names[name] = 1
            block_name = name
        elseif name != prev_name
            # New block os segments which has the same name as an older block.
            names[name] += 1
            block_name = name * string(names[name])
        end
        prev_name = name

        chromosome = parse(Int64, chromosome)
        startloc = parse(Int64, startloc)
        endloc = parse(Int64, endloc)
        push!(df, [block_name, chromosome, startloc, endloc])
    end
    return df
end

"""
    remove_false_matches(child_matches, parent_matches)

Remove matches from the child_matches if corresponding parent_matches
share less cM with the match than the child.

Return the matches where both parents have contributed and
the matches that are from one parent only.
"""
function remove_false_matches(child_matches::DataFrame, parent_matches::DataFrame)
    both = similar(child_matches, 0)
    parent_cleaned = similar(child_matches, 0)
    # Create a dictionary of parent matches cM.
    parentCM = Dict{String, Float64}()
    for row in eachrow(parent_matches)
        name = row[:Name]
        cM = row[:Total_cM_shared]
        parentCM[name] = cM
    end
    # From the matches that were inherited from the parent
    # remove thosw where child shared cM > parent shared cM.
    for row in eachrow(child_matches)
        name = row[:Name]
        childCM = row[:Total_cM_shared]
        if childCM > parentCM[name]
            push!(both, row)
        else
            push!(parent_cleaned, row)
        end
    end
    return both, parent_cleaned
end

"""
    most_common_ethnicity(d::Dict{String, Integer})::String

Eetermine the ethnicity based on a map that contains a count
of matches from each country.

d is a dictionary that contains the number of matches for each country.

Returns the name of the country that harbors the most matches
or an empty string if there is no clear result.
"""
function most_common_ethnicity(d::Dict{String, Integer})::String
    if length(keys(d)) == 0
        return ""
    end

    # Determine country with the most cousins.
    country = ""
    top = 0
    second = 0
    for (key, value) in d
        DEBUG && print(key, " ", value, "  ")
        if value > second && value <= top
            second = value
        end
        if value > top
            country = key
            second = top
            top = value
        end
    end
    DEBUG && println()

    # Look if the most common country is valid.
    # Experiments with different conditions are allowed here.

    # All matches from one country.
    if second == 0
        return country
    end

    # Condition for the top most matches to be valid.
    if top / second >= 1.5
        return country 
    else
        return ""
    end
end



