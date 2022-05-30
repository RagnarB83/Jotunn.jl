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
#@timeit to "jSCF(compile_4cloop)" result_4c=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="4c")
#for i in 1:repeat
#    @timeit to "jSCF(4cloop)" jSCF(H2O, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="4c")
#end

#4c-turbo
#@timeit to "jSCF(compile_4cturbo)" result_4ct=jSCF(mol, basis; maxiter=200, fock_algorithm="turbo", WFtype="RHF", tei_type="4c")
#for i in 1:repeat
#    @timeit to "jSCF(4cturbo)" jSCF(H2O, basis; maxiter=200, fock_algorithm="turbo", WFtype="RHF", tei_type="4c")
#end


#sp-loop bad
#@timeit to "jSCF(compile_sp_old)" result_sp_old=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")
#for i in 1:repeat
#    @timeit to "jSCF(sp1)" jSCF(H2O, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")
#end

#sp-loop good
@timeit to "jSCF(compile_sp)" result_sp=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")
for i in 1:repeat
    @timeit to "jSCF(sp)" jSCF(H2O, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")
end

#println("Energy (4c): $(result_4c["energy"])")
#println("Energy (4ct): $(result_4ct["energy"])")
#println("Energy (sp_old): $(result_sp_old["energy"])")
println("Energy (sp): $(result_sp["energy"])")


#Print timing output
show(to)


