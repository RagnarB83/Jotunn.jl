using Jotunn


Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)
#SCANNING VALUES
lshift_values=[0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,5.0,10.0,100.0]
damping_values=[0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
damping_thresh_values=[1000,100,10,1,1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7]
lshift_thresh_values=[1000,100,10,1,1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7]
damping_val=0.4


###
#TODO: Automate further once we have DIIS

#Testing different damping values
######################################
results=[]
fixed_damping_thresh=0.001
for VAL in damping_values
    result= jSCF(Acetic, "sto-3g"; WFtype="UHF", damping=true, levelshift=false, damping_val=VAL, damping_thresh=fixed_damping_thresh, maxiter=600)
    finaliter=result["finaliter"]
    finalenergy=result["energy"]
    push!(results,"Damping value: $VAL . Iterations: $finaliter Energy: $finalenergy")
end
println("FINAL DAMPING VALRESULTS\n")
println("Fixed fixed_damping_thresh: $fixed_damping_thresh")
for r in results
    println(r)
end
######################################

#Testing different damping thresholds
######################################
results=[]
fixed_damping_val=0.2
for VAL in damping_thresh_values
    result= jSCF(Acetic, "sto-3g"; WFtype="UHF", damping=true, levelshift=false, damping_val=fixed_damping_val, damping_thresh=VAL, maxiter=600)
    finaliter=result["finaliter"]
    finalenergy=result["energy"]
    push!(results,"Damping thresh: $VAL . Iterations: $finaliter Energy: $finalenergy")
end
println("FINAL DAMPING THRESH RESULTS\n")
println("Fixed damping val: $fixed_damping_val")
for r in results
    println(r)
end
######################################


#Testing different levelshift values
######################################
results=[]
fixed_lshift_thresh=0.001
for VAL in lshift_values
    result= jSCF(Acetic, "sto-3g"; WFtype="UHF", damping=false, levelshift=true, lshift_thresh=fixed_lshift_thresh, levelshift_val=VAL, maxiter=600)
    finaliter=result["finaliter"]
    finalenergy=result["energy"]
    push!(results,"Levelshift value: $VAL . Iterations: $finaliter Energy: $finalenergy")
end
println("FINAL Leveslhift VALUE RESULTS\n")
println("Fixed lshift thresh : $fixed_lshift_thresh")
for r in results
    println(r)
end
######################################

#Testing different levelshift thresh values
######################################
results=[]
fixed_levelshift_val=0.4
for VAL in lshift_thresh_values
    result= jSCF(Acetic, "sto-3g"; WFtype="UHF", damping=false, levelshift=true, levelshift_val=fixed_levelshift_val, lshift_thresh=VAL, maxiter=600)
    finaliter=result["finaliter"]
    finalenergy=result["energy"]
    push!(results,"Levelshift thresh: $VAL . Iterations: $finaliter Energy: $finalenergy")
end
println("FINAL Leveslhift RESULTS\n")
println("Fixed lshift val: $fixed_levelshift_val")
for r in results
    println(r)
end
######################################
