# Jotunn: a simple quantum chemistry program in Julia
RHF/UHF code written in in Julia (with integrals provided by libcint via GaussianBasis.jl). 
Primarily written to to be a simple, but general and easy-to-understand RHF/UHF code without being horribly slow.
Might one day turn into something useful.



## Current features
- 1 and 2-electron integrals via [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) and [libcint](https://github.com/sunqm/libcint).
- Conventional RHF and UHF algorithm with full rank-4 tensor.
- Fock matrix speedup via [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
- Mulliken population analysis (RHF and UHF)
- Mayer bond orders (RHF and UHF)
- SCF convergence aids: levelshifting
- Basis sets:
    - Support for all internal basis sets in [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl).
    - Support for reading in external basis set in ORCA format.

## Development: features to be added
- Improving speed: 
    - Support sparse 2-electron integrals
    - Density fitting
- Improving SCF convergence:
    - DIIS
    - Damping 
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

## Development: bugs to be fixed
- Fix bug: GaussianBasis.jl sensitive to scientific notation in coordinates
- Read-in basis set: Normalize contraction coefficients

## Dependencies (not in Julia standard library)
- [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) (1 and 2-electron integrals)
- [Molecules.jl](https://github.com/FermiQC/Molecules.jl) (helper program with GaussianBasis.jl)
- [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) (to speed up Fock matrix)
- [PrettyTables.jl](https://github.com/ronisbr/PrettyTables.jl) (pretty output)
- [Crayons.jl](https://github.com/KristofferC/Crayons.jl) (pretty output)


## Documentation:


### How to install:

*Option 1. Manual setup:*

First clone or download the Jotunn source-code. 

Make package available to Julia by setting the environment variable:
```sh
export JULIA_LOAD_PATH=/path/to/Jotunn/src:$JULIA_LOAD_PATH  
#copy-paste in Unix shell to make Jotunn package available to Julia
```

*Option 2. Install package:*

Launch a Julia session and copy-paste:
```julia
using Pkg; Pkg.add(url="https://github.com/RagnarB83/Jotunn.jl")
```

#
### How to use:

The primary Jotunn functions to use are **create_fragment** (a Jotunn molecule fragment object) and 
**jHF** (the Jotunn RHF/UHF code).

*Option 1: 
Launch an interactive julia session:*
```sh
julia
```

Then within Julia REPL, import Jotunn, create a molecule fragment and call jHF :
```julia
using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
jHF(H2, "sto-3g")
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

See  See [GaussianBasis/lib directory](https://github.com/FermiQC/GaussianBasis.jl/tree/main/lib) for list of available basis-sets.

**simple-input.jl:**
```julia
using Jotunn

# Create molecular fragments

H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#Simple call
energy= jHF(H2, "sto-3g")
println("Energy from Jotunn: $energy Eh")
 ```
**moreoptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More keywords
energy= jHF(H2O, "sto-3g"; maxiter=200, fock_algorithm="turbo", printlevel=2,
    HFtype="RHF", levelshift=2.0, lshift_thresh=1e-4, tei_type="4c", 
    print_final_matrices=true, debugprint=true)
 ```

**alloptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#All features
energy = jHF(H2O, "sto-3g"; HFtype="UHF", guess="hcore", 
    basisfile="none", maxiter=200, printlevel=1,
    print_final_matrices=false, debugprint=false, 
    rmsDP_threshold=5e-9, maxDP_threshold=1e-7, energythreshold=1e-8, 
    tei_type="4c", fock_algorithm="turbo", 
    levelshift=1.0, lshift_thresh=0.001)
 ```

#
### Example output:

```text

==================================================
                   JOTUNN
a simple quantum chemistry program in Julia
==================================================

jHF module: a RHF/UHF program

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
  1.763071 seconds (3.59 M allocations: 169.882 MiB, 3.40% gc time, 99.00% compilation time)
Choosing Fock algorithm.
Small system (Basis dim: 24). Choosing loop Fock.
Providing guess for density matrix

CALCULATION SETTINGS
====================== =================
              HF type              RHF  
            Basis set         def2-svp  
  No. basis functions               24  
                Guess            hcore  
      2-electron type               4c  
       Fock algorithm             loop  
  S lowest eigenvalue           0.0281  
====================== =================

Beginning SCF iterations
  Iter            Energy            deltaE            RMS-DP            Max-DP Levelshift    Damping
     1   -121.4976082145   -121.4976082145      0.3526149661      0.8363810730       true         no
     2    -67.0737370632     54.4238711514      0.2083579545      2.6274265828       true         no
     3    -74.8306068274     -7.7568697642      0.1104459981      1.2059710798       true         no
     4    -73.3621212492      1.4684855782      0.0145349590      0.1216172035       true         no
     5    -73.4905861685     -0.1284649193      0.0086857292      0.0960667848       true         no
     6    -73.4651054446      0.0254807239      0.0032598910      0.0345616540       true         no
     7    -73.4627145916      0.0023908530      0.0014452233      0.0154176001       true         no
     8    -73.4607371246      0.0019774670      0.0006174810      0.0065503619       true         no
     9    -73.4599903466      0.0007467780      0.0002679091      0.0028412302       true         no
    10    -73.4596529564      0.0003373902      0.0001159732      0.0012289148       true         no
    11    -73.4595082945      0.0001446619      0.0000502891      0.0005327990       true         no
    12    -72.1655274220      1.2939808725      0.2960680447      3.2213815245      false         no
    13    -78.6056576801     -6.4401302581      0.0970784541      0.5127348840      false         no
    14    -74.3290472853      4.2766103948      0.0518358027      0.6226871994      false         no
    15    -76.9591341296     -2.6300868443      0.0312625824      0.2042331625      false         no
    16    -75.3778166693      1.5813174603      0.0168315378      0.2065614086      false         no
    17    -76.3068587055     -0.9290420362      0.0101305710      0.0551039713      false         no
    18    -75.7595158677      0.5473428378      0.0056294333      0.0689660519      false         no
    19    -76.0788083030     -0.3192924353      0.0033640695      0.0185103155      false         no
    20    -75.8916698982      0.1871384048      0.0019082419      0.0231525213      false         no
    21    -76.0008069394     -0.1091370411      0.0011325964      0.0063937221      false         no
    22    -75.9369524458      0.0638544936      0.0006501422      0.0078102436      false         no
    23    -75.9742092353     -0.0372567895      0.0003838726      0.0021914433      false         no
    24    -75.9524269979      0.0217822374      0.0002218711      0.0026444978      false         no
    25    -75.9651410903     -0.0127140924      0.0001305227      0.0007492345      false         no
    26    -75.9577105969      0.0074304934      0.0000757456      0.0008977448      false         no
    27    -75.9620488171     -0.0043382202      0.0000444492      0.0002559043      false         no
    28    -75.9595139815      0.0025348356      0.0000258581      0.0003052966      false         no
    29    -75.9609941594     -0.0014801779      0.0000151492      0.0000873643      false         no
    30    -75.9601294002      0.0008647593      0.0000088262      0.0001039416      false         no
    31    -75.9606344143     -0.0005050141      0.0000051654      0.0000298180      false         no
    32    -75.9603393956      0.0002950186      0.0000030122      0.0000354144      false         no
    33    -75.9605116958     -0.0001723001      0.0000017616      0.0000101756      false         no
    34    -75.9604110468      0.0001006490      0.0000010279      0.0000120720      false         no
    35    -75.9604698313     -0.0000587846      0.0000006009      0.0000034721      false         no
    36    -75.9604354935      0.0000343378      0.0000003507      0.0000041164      false         no
    37    -75.9604555492     -0.0000200557      0.0000002050      0.0000011847      false         no
    38    -75.9604438343      0.0000117149      0.0000001197      0.0000014039      false         no
    39    -75.9604506768     -0.0000068424      0.0000000699      0.0000004042      false         no
    40    -75.9604466800      0.0000039967      0.0000000408      0.0000004789      false         no
    41    -75.9604490145     -0.0000023344      0.0000000239      0.0000001379      false         no
    42    -75.9604476509      0.0000013636      0.0000000139      0.0000001633      false         no
    43    -75.9604484473     -0.0000007964      0.0000000081      0.0000000471      false         no
    44    -75.9604479821      0.0000004652      0.0000000048      0.0000000557      false         no
    45    -75.9604482539     -0.0000002717      0.0000000028      0.0000000161      false         no
    46    -75.9604480951      0.0000001587      0.0000000016      0.0000000190      false         no
    47    -75.9604481878     -0.0000000927      0.0000000009      0.0000000055      false         no
    48    -75.9604481337      0.0000000541      0.0000000006      0.0000000065      false         no
    49    -75.9604481653     -0.0000000316      0.0000000003      0.0000000019      false         no
    50    -75.9604481469      0.0000000185      0.0000000002      0.0000000022      false         no
    51    -75.9604481576     -0.0000000108      0.0000000001      0.0000000006      false         no
    52    -75.9604481513      0.0000000063      0.0000000001      0.0000000008      false         no

                              SCF converged in 52 iterations! Hell yeah! 🎉

┌──────────────────────┬────────────────┐
│ Energy contributions │          E(Eh) │
├──────────────────────┼────────────────┤
│         Total energy │   -75.96044815 │
│    Nuclear repulsion │     9.11937991 │
│    Electronic energy │   -85.07982806 │
│    1-electron energy │  -122.90608314 │
│    2-electron energy │    37.82625509 │
│       Kinetic energy │    75.74591703 │
│     Potential energy │  -151.70636518 │
├──────────────────────┼────────────────┤
│         Virial ratio │     2.00283225 │
└──────────────────────┴────────────────┘
  1.078943 seconds (1.77 M allocations: 95.802 MiB, 2.10% gc time, 90.81% compilation time)


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
      Final HF energy     -75.96044815  
     Molecule formula              OHH  
  Number of electrons               10  
            Basis set         def2-svp  
              HF type              RHF  
       Fock algorithm             loop  
       SCF iterations               52  
====================== =================
  ```