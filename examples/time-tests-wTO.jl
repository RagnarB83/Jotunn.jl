using Jotunn
using TimerOutputs

#Simple time-tests using TimerOutputs

#Creating timer object
const to = TimerOutput()

#Timing create_fragment
@timeit to "create_fragment" H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)


basis="def2-qzvpp"

#Measure 1st run separately as it includes compilation
@timeit to "jHF(compile_turbo)" jHF(H2O, basis; maxiter=200, fock_algorithm="turbo", HFtype="RHF")
@timeit to "jHF(compile_loop)" jHF(H2O, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF")

#Run multiple times
for i in 1:3
    @timeit to "jHF(turbo)" jHF(H2O, basis; maxiter=200, fock_algorithm="turbo", HFtype="RHF")
end
#for i in 1:3
#    @timeit to "jHF(loop)" jHF(H2O, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF")
#end


#Print timing output
show(to)


