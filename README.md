# Jotunn: a simple quantum chemistry program in Julia
RHF/UHF code written in in Julia. 
Primarily written to to be a simple, but general and easy-to-understand RHF/UHF code without being horribly slow.
Might one day turn into something useful.



## Current features
- 1 and 2-electron integrals via [GaussianBasis.jl](https://github.com/FermiQC/GaussianBasis.jl) and [libcint](https://github.com/sunqm/libcint).
- RHF and UHF algorithm with full 4-rank tensor. 
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

*NOT YET READY...*

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
jHF(H2, "STO-3G")
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
energy= jHF(H2, "STO-3G")
println("Energy from Jotunn: $energy Eh")
 ```
**moreoptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More keywords
energy= jHF(H2O, "STO-3G"; maxiter=200, fock_algorithm="turbo", 
    HFtype="RHF", levelshift=2.0, lshift_thresh=1e-4, tei_type="4c", 
    print_final_matrices=true, debugprint=true)
 ```

**alloptions-input.jl:**
 ```julia
using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#All features
energy = jHF(H2O, "STO-3G"; HFtype="UHF", guess="hcore", 
    basisfile="none", maxiter=200, 
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
         Number of atoms                2  
        Molecule formula               HH  
                  Charge                0  
       Spin multiplicity                1  
           No. electrons                2  
  No. unpaired electrons                0  
       Nuclear repulsion       0.71510434  
========================= =================
Direct calculation of integrals via GaussianBasis.jl library
Creating basis set object
Looking up basis-name in GaussianBasis.jl

Calculating 1-electron integrals
Calculating 2-electron integrals
Choosing Fock algorithm.
Small system (Basis dim: 2). Choosing loop Fock.
Providing guess for density matrix

CALCULATION SETUP
====================== =================
              HF type              RHF  
            Basis set           STO-3G  
  No. basis functions       2.00000000  
                Guess            hcore  
             tei_type               4c  
       Fock algorithm             loop  
  S lowest eigenvalue       0.34012688  
====================== =================

Beginning SCF iterations

==============================
SCF iteration 1
==============================
Levelshift: true
Level shift active. Adding levelshift: 2.0 in iteration 1
Current energy: -1.586873896382793
Energy change: -1.586873896382793 Eh (threshold: 1.0e-8)
RMS-DP: 0.8825776863627137 (threshold: 5.0e-9)
Max-DP: -0.018125080890486533 (threshold: 1.0e-7)

==============================
SCF iteration 2
==============================
Levelshift: true
Level shift active. Adding levelshift: 2.0 in iteration 2
Current energy: -0.8491820078453758
Energy change: 0.7376918885374173 Eh (threshold: 1.0e-8)
RMS-DP: 0.04745591704180255 (threshold: 5.0e-9)
Max-DP: 0.07012163363267665 (threshold: 1.0e-7)

...

==============================
SCF iteration 14
==============================
Levelshift: false
Current energy: -1.1167593078220257
Energy change: 0.0 Eh (threshold: 1.0e-8)
RMS-DP: 4.1331225358956274e-9 (threshold: 5.0e-9)
Max-DP: 5.845118056235776e-9 (threshold: 1.0e-7)

SCF converged in 14 iterations! Hell yeah! 🎉
┌──────────────────────┬────────────────┐
│ Energy contributions │          E(Eh) │
├──────────────────────┼────────────────┤
│         Total energy │    -1.11675931 │
│    Nuclear repulsion │     0.71510434 │
│    Electronic energy │    -1.83186365 │
│    1-electron energy │    -2.50661957 │
│    2-electron energy │     0.67475593 │
│       Kinetic energy │     1.20128718 │
│     Potential energy │    -2.31804648 │
├──────────────────────┼────────────────┤
│         Virial ratio │     1.92963559 │
└──────────────────────┴────────────────┘

------------------------------
MO Energies (closed-shell)
------------------------------
                        ⍺
┌────┬──────┬──────────────┬──────────────┐
│ MO │ Occ. │        E(Eh) │        E(eV) │
├────┼──────┼──────────────┼──────────────┤
│  1 │  2.0 │      -0.5786 │     -15.7433 │
│  2 │  0.0 │       0.6711 │      18.2628 │
└────┴──────┴──────────────┴──────────────┘


Mulliken Population Analysis
┌──────┬─────────┬────────────┐
│ Atom │ Element │     Charge │
├──────┼─────────┼────────────┤
│    1 │       H │  -0.000000 │
│    2 │       H │   0.000000 │
└──────┴─────────┴────────────┘

Sum of charges: 0.0000

Mayer bond orders (>0.01):
┌─────────┬────────┐
│    Bond │    MBO │
├─────────┼────────┤
│ H1 - H2 │ 1.0000 │
└─────────┴────────┘

FINAL RESULTS
====================== =================
      Final HF energy      -1.11675931  
     Molecule formula               HH  
  Number of electrons       2.00000000  
            Basis set           STO-3G  
              HF type              RHF  
       Fock algorithm             loop  
       SCF iterations               14  
====================== =================
  ```