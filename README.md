# Jotunn: a simple quantum chemistry program in Julia
RHF/UHF code written in in Julia (with integrals provided by libcint via GaussianBasis.jl). 
Primarily written to to be a simple, but general and easy-to-understand RHF/UHF code without being horribly slow.
Might one day turn into something useful.


## Current features
- 1 and 2-electron integrals via [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) and [libcint](https://github.com/sunqm/libcint).
- Conventional RHF and UHF algorithm using 2-el integrals as full rank-4 tensor (4c) or sparse version (sparse4c)
- Fock matrix speedup (4c) via [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
- Mulliken population analysis (RHF and UHF)
- Mayer bond orders (RHF and UHF)
- SCF convergence aids: DIIS, levelshifting, static damping
- Basis sets:
    - Support for all internal basis sets in [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl).
    - Support for reading in external basis set in ORCA format.

## Development: features to be added
- Improving speed: 
    - Further fine-tuning of Fock code for sparse 2el-integral version.
    - Density fitting
- Improving SCF convergence:
    - Dynamic damping 
    - Better guess than Hcore
- DFT support:
    - interface to LibXC
    - DFT grids
    - hybrid-DFT
- Support broken-symmetry guess
- Simple electric properties (dipole, EFG)
- Hirshfeld population analysis
- Direct SCF algorithm ?
- Noncollinear HF/DFT
- Read-in basis set: Normalize contraction coefficients. Code currently assumes normalization.

## Dependencies (not in Julia standard library)
- [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) (1 and 2-electron integrals)
- [Molecules.jl](https://github.com/FermiQC/Molecules.jl) (helper program with GaussianBasis.jl)
- [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) (only to speed up Fock when using full 4-rank TEIs )
- [PrettyTables.jl](https://github.com/ronisbr/PrettyTables.jl) (pretty output)
- [Crayons.jl](https://github.com/KristofferC/Crayons.jl) (pretty output)


## Documentation:


### How to install:

*Option 1. Manual setup:*

First clone or download the Jotunn source-code. Make package available to Julia by setting the environment variable:
```sh
export JULIA_LOAD_PATH=/path/to/Jotunn/src:$JULIA_LOAD_PATH  
#copy-paste in Unix shell to make Jotunn package available to Julia
```

*Option 2. Install package:*

Launch a Julia session and copy-paste the line below into Julia REPL. This will install Jotunn as a Julia package:
```julia
using Pkg; Pkg.add(url="https://github.com/RagnarB83/Jotunn.jl")
```

#
### Basic functionality:



 **create_fragment** (a function to create a Jotunn molecule fragment object)
```julia
function create_fragment(;coords_string=nothing,xyzfile=nothing,pdbfile=nothing,fragfile=nothing, coords=nothing,
    elems=nothing, calc_connectivity=false, label=nothing, charge=nothing, mult=nothing)
```

**jSCF** (a function to run the Jotunn RHF/UHF code).

```julia
function jSCF(fragment, basisset="sto-3g"; WFtype::String="RHF", guess::String="hcore", basisfile::String="none", maxiter::Int64=120, 
    print_final_matrices::Bool=false, rmsDP_threshold::Float64=5e-9, maxDP_threshold::Float64=1e-7, tei_type::String="sparse4c",
    energythreshold::Float64=1e-8, debugprint::Bool=false, fock_algorithm::String="loop", 
    levelshift::Bool=false, levelshift_val::Float64=0.10, lshift_thresh::Float64=0.01,
    damping::Bool=true, damping_val::Float64=0.4, damping_thresh::Float64=0.01,
    diis::Bool=false, diis_size::Int64=7, diis_thresh::Float64=0.01,
    printlevel::Int64=1)
```

### How to use:

*Option 1: 
Launch an interactive julia session:*
```sh
julia
```

Then within Julia REPL, import Jotunn, create a molecule fragment and call jSCF :
```julia
using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
jSCF(H2, "sto-3g")
```


*Option 2: Create a Julia script (e.g. test.jl)*  



test.jl:
```julia
using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
```
Run inputscript like this: 
```sh
julia test.jl
```
#
### Example inputfiles:
*Example inputfiles below can all be found in examples directory*

See [GaussianBasis/lib directory](https://github.com/FermiQC/GaussianBasis.jl/tree/main/lib) for list of available basis-sets.

**simple-input.jl:**
```julia
using Jotunn

# Create molecular fragments

H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#Simple call
result= jSCF(H2, "sto-3g")
println("Result dictionary from Jotunn: $result")
println("Energy: $(result["energy"]) Eh")
 ```
**moreoptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More keywords
result= jSCF(H2O, "sto-3g"; maxiter=200, fock_algorithm="turbo", printlevel=2,
    WFtype="RHF", levelshift=2.0, lshift_thresh=1e-4, tei_type="4c", 
    print_final_matrices=true, debugprint=true)
 ```

**alloptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#All features
result = jSCF(H2O, basisset="STO-3G"; WFtype="RHF", guess="hcore", basisfile="none", maxiter=120, 
    print_final_matrices=false, rmsDP_threshold=5e-9, maxDP_threshold=1e-7, tei_type="4c",
    energythreshold=1e-8, debugprint=false, fock_algorithm="turbo", 
    levelshift=false, levelshift_val=0.10, lshift_thresh=0.01,
    damping=true, damping_val=0.4, damping_thresh=0.01,
    printlevel=1)
 ```

#
### Example output:

```text

==================================================
                   JOTUNN
a simple quantum chemistry program in Julia
==================================================

jSCF module: a RHF/UHF program

SYSTEM
========================= =================
         Number of atoms                3  
        Molecule formula              OHH  
                  Charge                0  
       Spin multiplicity                1  
           No. electrons               10  
  No. unpaired electrons                0  
       Nuclear repulsion       9.11937991  
========================= =================
Integrals provided via GaussianBasis.jl library
Creating basis set object

Calculating 1-electron integrals
Calculating 2-electron integrals
  2.461401 seconds (3.59 M allocations: 169.377 MiB, 3.08% gc time, 98.97% compilation time)
Choosing Fock algorithm.
Small system (Basis dim: 24). Choosing loop Fock.
Providing guess for density matrix
Energy of guess: 9.119379905339786 Eh

CALCULATION SETTINGS
========================= =================
                 HF type              RHF  
               Basis set         def2-svp  
     No. basis functions               24  
                   Guess            hcore  
         2-electron type               4c  
          Fock algorithm             loop  
     S lowest eigenvalue           0.0281  
              Levelshift            false  
    Levelshift parameter           0.1000  
  Lshift-turnoff thresh.           0.0100  
                 Damping             true  
       Damping parameter           0.4000  
    Damp-turnoff thresh.           0.0100  
                    DIIS            false  
           DIIS vec size                7  
            DIIS thresh.           0.0100  
========================= =================

Beginning SCF iterations
Iter         Energy          deltaE        RMS-DP        Max-DP  Lshift    Damp    DIIS
   1 -126.165023244  -126.165023244   0.422360951   0.766341001   false    true   false
   2  -90.891462805    35.273560439   0.144192193   0.976274234   false    true   false
   3  -77.074577330    13.816885475   0.177725075   1.918743295   false    true   false
   4  -76.127629367     0.946947963   0.073956975   0.835446553   false    true   false
   5  -76.010001447     0.117627920   0.029057353   0.340558595   false    true   false
   6  -75.972851329     0.037150118   0.011633177   0.141278309   false    true   false
   7  -75.962058521     0.010792808   0.004814077   0.060224092   false    true   false
   8  -75.957733508     0.004325013   0.000795942   0.004387088   false   false   false
   9  -75.959936967    -0.002203459   0.000275068   0.002149001   false   false   false
  10  -75.960270463    -0.000333496   0.000103929   0.000945157   false   false   false
  11  -75.960395859    -0.000125396   0.000040334   0.000397913   false   false   false
  12  -75.960425046    -0.000029186   0.000015887   0.000166315   false   false   false
  13  -75.960441935    -0.000016889   0.000006333   0.000067903   false   false   false
  14  -75.960444414    -0.000002479   0.000002553   0.000028294   false   false   false
  15  -75.960447558    -0.000003143   0.000001038   0.000011428   false   false   false
  16  -75.960447455     0.000000103   0.000000427   0.000004844   false   false   false
  17  -75.960448170    -0.000000715   0.000000177   0.000001934   false   false   false
  18  -75.960447999     0.000000171   0.000000075   0.000000849   false   false   false
  19  -75.960448191    -0.000000192   0.000000032   0.000000330   false   false   false
  20  -75.960448113     0.000000078   0.000000014   0.000000154   false   false   false
  21  -75.960448171    -0.000000058   0.000000006   0.000000057   false   false   false
  22  -75.960448141     0.000000030   0.000000003   0.000000029   false   false   false
  23  -75.960448160    -0.000000019   0.000000001   0.000000010   false   false   false
  24  -75.960448150     0.000000010   0.000000001   0.000000006   false   false   false
  25  -75.960448156    -0.000000006   0.000000000   0.000000002   false   false   false

                              SCF converged in 25 iterations! Hell yeah! 🎉

┌──────────────────────┬────────────────┐
│ Energy contributions │          E(Eh) │
├──────────────────────┼────────────────┤
│         Total energy │   -75.96044816 │
│    Nuclear repulsion │     9.11937991 │
│    Electronic energy │   -85.07982806 │
│    1-electron energy │  -122.90608315 │
│    2-electron energy │    37.82625508 │
│       Kinetic energy │    75.74591703 │
│     Potential energy │  -151.70636519 │
├──────────────────────┼────────────────┤
│         Virial ratio │     2.00283225 │
└──────────────────────┴────────────────┘
  1.332684 seconds (1.72 M allocations: 90.361 MiB, 2.09% gc time, 94.01% compilation time)


MO Energies (closed-shell)

┌────┬──────┬──────────────┬──────────────┐
│ MO │ Occ. │        E(Eh) │        E(eV) │
├────┼──────┼──────────────┼──────────────┤
│  1 │  2.0 │     -20.5475 │    -559.1252 │
│  2 │  2.0 │      -1.3154 │     -35.7944 │
│  3 │  2.0 │      -0.6979 │     -18.9900 │
│  4 │  2.0 │      -0.5683 │     -15.4630 │
│  5 │  2.0 │      -0.4978 │     -13.5452 │
│  6 │  0.0 │       0.1748 │       4.7574 │
│  7 │  0.0 │       0.2542 │       6.9167 │
│  8 │  0.0 │       0.7862 │      21.3927 │
│  9 │  0.0 │       0.8642 │      23.5149 │
│ 10 │  0.0 │       1.1852 │      32.2514 │
│ 11 │  0.0 │       1.2017 │      32.7013 │
│ 12 │  0.0 │       1.2705 │      34.5727 │
│ 13 │  0.0 │       1.3419 │      36.5163 │
│ 14 │  0.0 │       1.5981 │      43.4873 │
│ 15 │  0.0 │       1.6553 │      45.0430 │
│ 16 │  0.0 │       1.8029 │      49.0589 │
│ 17 │  0.0 │       2.0563 │      55.9553 │
│ 18 │  0.0 │       2.5418 │      69.1649 │
│ 19 │  0.0 │       2.5863 │      70.3762 │
│ 20 │  0.0 │       3.3239 │      90.4486 │
│ 21 │  0.0 │       3.3739 │      91.8083 │
│ 22 │  0.0 │       3.5635 │      96.9678 │
│ 23 │  0.0 │       3.8870 │     105.7716 │
│ 24 │  0.0 │       4.2295 │     115.0902 │
└────┴──────┴──────────────┴──────────────┘

Mulliken Population Analysis (closed-shell)
┌──────┬─────────┬────────────┐
│ Atom │ Element │     Charge │
├──────┼─────────┼────────────┤
│    1 │       O │  -0.353914 │
│    2 │       H │   0.176957 │
│    3 │       H │   0.176957 │
└──────┴─────────┴────────────┘

Sum of charges: 0.0000

Mayer bond orders
Threshold: 0.01
┌─────────┬────────┐
│    Bond │    MBO │
├─────────┼────────┤
│ O1 - H2 │ 0.9947 │
│ O1 - H3 │ 0.9947 │
└─────────┴────────┘

FINAL RESULTS
====================== =================
      Final HF energy     -75.96044816  
     Molecule formula              OHH  
  Number of electrons               10  
            Basis set         def2-svp  
              HF type              RHF  
       Fock algorithm             loop  
       SCF iterations               25  
====================== =================
  ```
