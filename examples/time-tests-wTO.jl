using Jotunn
using TimerOutputs

#Simple time-tests using TimerOutputs

#Creating timer object
const to = TimerOutput()

#Timing create_fragment
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
mol=H2O

basis="def2-tzvpp"

repeat=5

#4c-loop
#@timeit to "jHF(compile_4cloop)" result_4c=jHF(mol, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="4c")
#for i in 1:repeat
#    @timeit to "jHF(4cloop)" jHF(H2O, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="4c")
#end

#4c-turbo
#@timeit to "jHF(compile_4cturbo)" result_4ct=jHF(mol, basis; maxiter=200, fock_algorithm="turbo", HFtype="RHF", tei_type="4c")
#for i in 1:repeat
#    @timeit to "jHF(4cturbo)" jHF(H2O, basis; maxiter=200, fock_algorithm="turbo", HFtype="RHF", tei_type="4c")
#end


#sp-loop bad
#@timeit to "jHF(compile_sp_old)" result_sp_old=jHF(mol, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")
#for i in 1:repeat
#    @timeit to "jHF(sp1)" jHF(H2O, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")
#end

#sp-loop good
@timeit to "jHF(compile_sp)" result_sp=jHF(mol, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")
for i in 1:repeat
    @timeit to "jHF(sp)" jHF(H2O, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")
end

#println("Energy (4c): $(result_4c["energy"])")
#println("Energy (4ct): $(result_4ct["energy"])")
#println("Energy (sp_old): $(result_sp_old["energy"])")
println("Energy (sp): $(result_sp["energy"])")


#Print timing output
show(to)


