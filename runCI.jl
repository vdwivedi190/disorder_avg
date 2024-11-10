using Dates, DelimitedFiles, Distributed, LinearAlgebra

# Run this as julia [PROJ_DIR] [CODE_DIR] from the terminal. 
# The latter parameter is optional. If it is missing, then the current directory is taken as the location of the julia code. 
# The PROJ_DIR must contain a job_list.txt generated by genscriptCI.jl

# Include the relevant Julia files
global codedir = length(ARGS) < 2 ? "./" : ARGS[2]
include( string(codedir,"/computeLyapunov.jl") )
include( string(codedir,"/computeCI.jl") ) 
println("Imported code.")


# Import the list of jobs from the dir given as the first argument. 
jobdir = ARGS[1]
try 
    cd(jobdir)
catch 
    println("The directory $jobdir does not exist!")
    exit() 
end 

jobfile = "job_list.txt"
try
	open(jobfile, "r")
catch 
	println("No job_list.txt file found in $jobdir. Aborting...")
	exit()
end 
global joblist = readdlm(jobfile,'\t',skipstart=1) 
global njobs = size(joblist)[1]
println("Found $njobs jobs in the file $jobdir/$jobfile")

# Store up the start time 
starting_at=now()

# Create the log file 
logfile = "run_log.txt" 
open(logfile, "w") do f
	 write(f,"Created the log file at $starting_at \n\n")
	 write(f,"jobID\tm\tW\tLy\t Nx\tworkerID\thostID\t\tTime taken\n")
	 println("Created the log file $jobdir/$logfile")
end

global qrflag = 0  # flag whether the results of QR decomposition are written to file
global uc = 2   # Dofs per unit cell = Number of bands
global dtype = 1  # Disorder type 

println("Computing disorder-averaged Lyapunov exponents for Chern insulator...")
println("Note that the additional setup means that first job might take unusually long.")

for job_index in 1:njobs
    println("\nWorking on input # ",job_index)
    perform_job(joblist[job_index,:], dtype, qrflag)
	flush(stdout)
end


stopping_at = now()
time_diff= Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(stopping_at) - Dates.DateTime(starting_at)))
println("Done. Time taken = $(time_diff)")