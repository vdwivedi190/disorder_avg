using Distributed, Dates, LinearAlgebra, DelimitedFiles

# Number of threads for basic linear algebra functions
BLAS.set_num_threads(4)

include("computeLyapunov.jl")
include("computeSOSL.jl")
jobfile = "./job_list.txt"   # Read the list of jobs from this file

starting_at=now()

println("Computing disorder-averaged Lyapunov exponents for the Kitaev spin liquid on the Shastry-Sutherland lattice with flux disorder.")


# # Parse the command string: julia sosl.jl ...
#  ==========================
#  1) Index of the task in the job file
#  2) Flag to write the Q and R matrices (1 for yes, default = no)
job_index =  parse(Int, ARGS[1])       # Index of the task
qrflag = length(ARGS) < 2 ? 0 : parse(Int, ARGS[2])

@everywhere uc = 4   # Degrees of freedom per unit cell = Number of bands

println("Working on input # ",job_index," in ",jobfile)
my_job = readdlm(jobfile,'\t',skipstart=1)[job_index,:]
perform_job(my_job,qrflag)
flush(stdout)

stopping_at=now()

time_diff= Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(stopping_at) - Dates.DateTime(starting_at)))
println("Done. Time taken = $(time_diff)")
