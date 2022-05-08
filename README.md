# Jotunn: a simple quantum chemistry program in Julia
RHF/UHF code written in in Julia. 
Primarily written to to be a simple, but general and easy-to-understand RHF/UHF code without being horribly slow.
Might one day turn into something useful.



**Current features**
- 1 and 2-electron integrals via GaussianBasis.jl and Libcint.
- RHF and UHF algorithm with full 4-rank tensor. 
- Mulliken population analysis (RHF and UHF)
- Mayer bond orders (RHF and UHF)
- SCF convergence aids: levelshifting
- Basis sets:
    - Support for all internal basis sets in GaussianBasis.
    - Support for reading in external basis set in ORCA format.

**Development: features to be added**
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

**Development: bugs to be fixed **
- Fix bug: GaussianBasis.jl sensitive to scientific notation in coordinates
- Read-in basis set: Normalize contraction coefficients


**Documentation:**


**Example:**

#
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