import sys
import numpy as np

def symmetrise(Mat): # Symmetrize a matrix given a triangular one
    return Mat + Mat.T - np.diag(Mat.diagonal())

def eint(a,b,c,d): # Return compound index given four indices using Yoshimine sort
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab

    return abcd

def tei(a, b, c, d): # Return value of two electron integral
    return twoe.get(eint(a, b, c, d), 0) 

def fprime(X, F): # Put Fock matrix in orthonormal AO basis
    return np.dot(np.transpose(X), np.dot(F, X)) 

def makedensity(C, D, dim, Nelec): # Make density matrix and store old one to test for convergence
    Dold = np.zeros((dim, dim))
    for mu in range(0, dim):
        for nu in range(0, dim):
            Dold[mu,nu] = D[mu, nu]
            D[mu,nu] = 0
            for m in range(0, int(Nelec/2)):
                D[mu,nu] = D[mu,nu] + 2*C[mu,m]*C[nu,m]

    return D, Dold 

def makefock(Hcore, P, dim): # Make Fock Matrix
    F = np.zeros((dim, dim))
    for i in range(0, dim):
        for j in range(0, dim):
            F[i,j] = Hcore[i,j]
            for k in range(0, dim):
                for l in range(0, dim):
                    F[i,j] = F[i,j] + P[k,l]*(tei(i+1,j+1,k+1,l+1)-0.5*tei(i+1,k+1,j+1,l+1))
    
    return F 

def deltap(D, Dold): # Calculate change in density matrix using Root Mean Square Deviation (RMSD)
    DELTA = 0.0
    for i in range(0, dim):
        for j in range(0, dim):
            DELTA = DELTA + ((D[i,j] - Dold[i,j])**2)

    return (DELTA)**(0.5)

def currentenergy(D, Hcore, F, dim): # Calculate energy at iteration
    EN = 0
    for mu in range(0, dim):
        for nu in range(0, dim):
            EN += 0.5*D[mu,nu]*(Hcore[mu,nu] + F[mu,nu])
            
    return EN

Nelec = 2 # The number of electrons in our system 
ENUC = np.genfromtxt('./enuc.dat',dtype=float, delimiter=',') # ENUC = nuclear repulsion, 
Sraw = np.genfromtxt('./sint.dat',dtype=None)                    # Sraw is overlap matrix, 
Traw = np.genfromtxt('./tint.dat',dtype=None)                    # Traw is kinetic energy matrix,
Vraw = np.genfromtxt('./vint.dat',dtype=None)                    # Vraw is potential energy matrix

dim = 2 # dim is the number of basis functions 
S = np.zeros((dim, dim)) # Initialize integrals, and put them in a Numpy array
T = np.zeros((dim, dim))
V = np.zeros((dim, dim))

for i in Sraw: S[i[0]-1, i[1]-1] = i[2] # Put the integrals into a matrix 
for i in Traw: T[i[0]-1, i[1]-1] = i[2] # Put the integrals into a matrix
for i in Vraw: V[i[0]-1, i[1]-1] = i[2] # Put the integrals into a matrix

S            = symmetrise(S) # Flip the triangular matrix in the diagonal 
V            = symmetrise(V) # Flip the triangular matrix in the diagonal
T            = symmetrise(T) # Flip the triangular matrix in the diagonal
TEI          = np.genfromtxt('./two_elec_int.dat') # Load two electron integrals
twoe         = {eint(row[0], row[1], row[2], row[3]) : row[4] for row in TEI} # Put in python dictionary
Hcore        = T + V # Form core Hamiltonian matrix as sum of one electron kinetic energy, T and potential energy, V matrices
SVAL, SVEC   = np.linalg.eigh(S) # Diagonalize basis using symmetric orthogonalization 
SVAL_minhalf = (np.diag(SVAL**(-0.5))) # Inverse square root of eigenvalues
S_minhalf    = np.dot(SVEC, np.dot(SVAL_minhalf, np.transpose(SVEC)))
P            = np.zeros((dim, dim)) # P represents the density matrix, Initially set to zero.
DELTA        = 1 # Set placeholder value for delta
count        = 0 # Count how many SCF cycles are done, N(SCF)

while DELTA > 0.00001:
    count     += 1                             # Add one to number of SCF cycles counter
    F         = makefock(Hcore, P, dim)        # Calculate Fock matrix, F
    Fprime    = fprime(S_minhalf, F)           # Calculate transformed Fock matrix, F'
    E, Cprime = np.linalg.eigh(Fprime)         # Diagonalize F' matrix
    C         = np.dot(S_minhalf, Cprime)      # 'Back transform' the coefficients into original basis using transformation matrix
    P, OLDP   = makedensity(C, P, dim, Nelec)  # Make density matrix
    DELTA     = deltap(P, OLDP)                # Test for convergence. If criteria is met exit loop and calculate properties of interest
    
    print("E = {:.6f}, N(SCF) = {}".format(currentenergy(P, Hcore, F, dim) + ENUC, count))

print("SCF procedure complete, TOTAL E(SCF) = {} hartrees".format(currentenergy(P, Hcore, F, dim) + ENUC))