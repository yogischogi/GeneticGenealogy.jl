# GeneticGenealogy

GeneticGenealogy is a Julia module that contains methods
to

1. calculate ethnicity estimates from MyHeritage matches.
2. calculate TMRCAs on a phylogenetic tree.


## Installation and usage

1. Install [Julia](https://julialang.org/downloads/) on
   your computer.
2. At the Julia prompt type `]` to switch to the package
   manager. Type `add https://github.com/yogischogi/GeneticGenealogy.jl`
   to install the package.
3. Type `using GeneticGenealogy` to use the package.
4. You can now use the `ethnicity` function to calculate ethnicity
   estimates.

## Calculate ethnicity estimate

To run this example you need your `matches.csv` and your `shared_segments.csv`
file from MyHeritage. You can find them under `MyHeritage -> DNA Matches ->
The little menu with the dots at the right`. Important: You need the English
versions!

```julia
julia> ]
pkg> add https://github.com/yogischogi/GeneticGenealogy.jl
julia> using GeneticGenealogy
julia> ethnicity("matches.csv", "shared_segments.csv", birth_country="Germany")
Total: 1124 segments
        Germany  710  63.2%
        Denmark  164  14.6%
    Netherlands  106   9.4%
```

### How does it work?

Most big companies calculate ethnicity estimates based on
reference populations. The result shows the genetic similarity to
their reference populations.

Very often this is not what the user expects. If you compare 
a North German to a Bavarian you will immediately realize that
they are different. There is no doubt that both of them are
German but not necessarily typical German and this is where the
trouble starts. How do you define a typical German? A North
German or a Bavarian may not think of themselves as typical German
and they may even be proud of it.

As a result companies often define a typical German as a German
who has SNP mutations that are common in Germany and less
common in other countries. But this also means that not every
German has them and indeed some Germans end up with having no
German ethnicity at all. No wonder that some bloggers claim
that ethnicity estimates are only good for party talk.

This package takes a different approach. It does not try to
calculate a genetic similarity. Instead it looks at the DNA
segments of your matches and determines in which country they
are most common. This works reasonable well unless a country
has experienced a mass migration in recent times. Such countries
can be excluded from the calculation.

The approach has the benefit that it does not need any reference
populations. If your ancestors were from England they were English,
regardless of whether they were typical English or not.

The downside of this method is that it needs many relatives for
a good ethnicity estimate. Also if some countries test more than
others or are not represented at all this distorts the results.


## Calculate TMRCAs on a phylogenetic tree

To calculate TMRCA estimates it is necessary to install the
[Phyloage](https://github.com/yogischogi/phyloage) program first.
Make sure that it is within your `PATH` variable so that it can be found.

You will need a file that contains a phylogenetic tree with
samples and the number of SNP mutations for each sample and each
branch. A simple tree looks like this:

```
P310
    P312, snps: 2
        id:001, snps: 30, Celtic warrior
        id:002, snps: 36, Spanish bull fighter
    U106, snps:1
        id:0815, snps: 34, Germanic hog hunter
        id:003, snps: 36, German waggon builder
```

Let the file be called `P310.txt`. Now you can calculate the
TMRCAs:

```julia
julia> using GeneticGenealogy
julia> a = readtree("P310.txt")
julia> b = calculateTMRCAs(a)
julia> writetree("results.txt", b)
```

The file `results.txt` contains the tree with TMRCA estimates
and confidence intervals.

```
P310, TMRCA: 5031, CI:[4249, 5959]
    P312, TMRCA: 4680, CI:[3674, 5966]
        id:001, Celtic warrior
        id:002, Spanish bull fighter
    U106, TMRCA: 4960, CI:[3920, 6280]
        id:0815, Germanic boar hunter
        id:003, German waggon builder
```

You can calibrate your tree by using different mutation rates.
In this case `calculateTMRCAs` allows an optional calibration parameter,
for example `calculateTMRCAs(a, cal=140.0)` if one mutation occurs
in 140 years on average.

It is also possible to merge several trees. This is very useful
if trees get very big and different persons work on different branches.

Imagine that the branches of our previous tree are stored in three
different files.

```
backbone.txt:
P310
    P312
    U106

P312.txt:
P312, snps: 2
    id:001, snps: 30, Celtic warrior
    id:002, snps: 36, Spanish bull fighter

U106.txt:
U106, snps:1
    id:0815, snps: 34, Germanic hog hunter
    id:003, snps: 36, German waggon builder
```

Merging them together is done by

```julia
julia> d = readtrees(["backbone.txt", "P312.txt", "U106.txt"])
```

The result is of course:

```
P310
    P312, snps: 2
        id:001, snps: 30, Celtic warrior
        id:002, snps: 36, Spanish bull fighter
    U106, snps:1
        id:0815, snps: 34, Germanic hog hunter
        id:003, snps: 36, German waggon builder
```






















