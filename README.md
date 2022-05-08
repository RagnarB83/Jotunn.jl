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

### Example:

#
**Julia inputfile.jl:**
```sh
include("../jHF.jl")

\# Create molecular fragments

H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

\#Simple call
energy= jHF(H2, "STO-3G")

\#More keywords
energy= jHF(H2, "STO-3G"; maxiter=200, fock_algorithm="turbo", HFtype="RHF", levelshift=2.0, lshift_thresh=1e-4, tei_type="4c", print_final_matrices=true, debugprint=true)

\#All features
energy = jHF(H2O, basisset="STO-3G"; HFtype="UHF", guess="hcore", basisfile="none", maxiter=200, 
    print_final_matrices=false, debugprint=false, 
    rmsDP_threshold=5e-9, maxDP_threshold=1e-7, energythreshold=1e-8, 
    tei_type="4c", fock_algorithm="turbo", 
    levelshift=1.0, lshift_thresh=0.001)


 ```

 **Output:**

```sh

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
  2.041039 seconds (3.49 M allocations: 155.288 MiB, 4.26% gc time, 99.82% compilation time)
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
  0.820469 seconds (776.32 k allocations: 41.955 MiB, 3.41% gc time, 76.26% compilation time)
  0.016378 seconds (14.32 k allocations: 807.283 KiB, 99.74% compilation time)
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
  0.011117 seconds (9.48 k allocations: 552.631 KiB, 98.86% compilation time)


Mulliken Population Analysis
┌──────┬─────────┬────────────┐
│ Atom │ Element │     Charge │
├──────┼─────────┼────────────┤
│    1 │       H │  -0.000000 │
│    2 │       H │   0.000000 │
└──────┴─────────┴────────────┘

Sum of charges: 0.0000
  0.154539 seconds (265.43 k allocations: 14.903 MiB, 98.16% compilation time)
  0.117040 seconds (7.99 k allocations: 300.141 KiB, 99.94% compilation time)

Mayer bond orders (>0.01):
┌─────────┬────────┐
│    Bond │    MBO │
├─────────┼────────┤
│ H1 - H2 │ 1.0000 │
└─────────┴────────┘
  0.119578 seconds (123.70 k allocations: 6.427 MiB, 99.02% compilation time)

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