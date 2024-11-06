using Distributed, Dates, LinearAlgebra, DelimitedFiles

# Number of threads for basic linear algebra functions
BLAS.set_num_threads(4)

include("computeLyapunov.jl")
include("computeCI.jl")
jobfile = "./job_list.txt"   # Read the list of jobs from this file

starting_at=now()
println("Computing disorder-averaged Lyapunov exponents for CHERN INSULATOR with on-site disorder.")

# Parse the command string: julia ci.jl ...
#  ==========================
#  1) Index of the task in the job file
#  2) Type of disorder: (default = 1)
#     1,2 = diagonal,sigma_z disorder on every site
#     3,4 = diagonal,sigma_z disorder with density rho, i.e, on rho*Ly unit cells chosen randomly for each iteration
#  3) Flag to write the Q and R matrices (1 for yes, default = no)
job_index =  parse(Int, ARGS[1])
dtype = length(ARGS) < 2 ? 1 : parse(Int, ARGS[2])
qrflag = length(ARGS) < 3 ? 0 : parse(Int, ARGS[2])

@everywhere uc = 2   # Degrees of freedom per unit cell = Number of bands

println("Working on input # ",job_index," in ",jobfile)
my_job = readdlm(jobfile,'\t',skipstart=1)[job_index,:]
perform_job(my_job,dtype,qrflag)
flush(stdout)

stopping_at=now()

time_diff= Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(stopping_at) - Dates.DateTime(starting_at)))
println("Done. Time taken = $(time_diff)")
