using Jotunn
#using BenchmarkTools


#Simple time-tests using elapsed
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

times=[]
#More options and more printing
for i in 1:5
    t = @elapsed jSCF(H2O, "sto-3g"; maxiter=200, fock_algorithm="loop", WFtype="RHF")
    push!(times,t)
end

println("times: $times")
