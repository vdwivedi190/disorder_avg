# This script generates the list of parameter values for which the disorder-averaging is to be performed and 
# puts it in a directory created for the project, along with subdirectories for run logs. It also optionally generates the bash scripts (one for each Ly) required to run it at a computing cluster. It was written for the CHEOPS cluster at the University of Cologne. 

# PARAMETERS
# ============
# The project name is also the name of the subdirectory that will be created. 
# The relevant system parameters here are J0, Jdelta, Jz 
# (See Physical Review Research 2 (4), 043159 (2020) for the definition of the first three parameters)
# The "disorder strength" (W) refers to the probability that a vison exists at each unit cell
# One can also scan over a range of the width (Ly) and the aspect ratio (scl = Lx/Ly)
# The variable pbc (= 0 or 1) sets the periodic boundary condition

projname = "test_sosl"

Lylist = [20 28 40 56 80]
sclist = [1]

J0list = 0.1:0.1:0.8
Jdeltalist = [0.0]
Jzlist = 1.5:0.1:1.8
Wlist = [0.50]   # Probability of a vison in each 4-plaquette
Elist = [0.0]
pbclist = [1]


bashflag = false  # Generate bash scripts for HPC cluster if true 


# =================================================================

# The following are relevant only for the bash script
ncpus = 4   
codepath = "./runSOSL_cluster.jl"

# Estimate the time required and format to a string suitable for the CHEOPS (Cologne HPC cluster) bash script
function computetime(Ly, scale)
    runtime = 1.85e-5 * Ly^2.9 * scale  # In seconds
	minutes = ceil(runtime/60)
	minutes = 10*Int64.(ceil(minutes/10))

    if minutes >= 60
        hours = Int64.(floor(minutes / 60))
        minutes = Int64.(minutes % 60)
	else
		hours = 0
    end

    return string("$hours:",string(minutes, pad=2),":00")
end


# Function to write out bash scripts for CHEOPS
function genscript(f,istart,iend,walltime)
    println(f,"#!/bin/bash -l\n")
    println(f,"#SBATCH --ntasks=1")
    println(f,"#SBATCH --cpus-per-task=$ncpus")
    println(f,"#SBATCH --array=",istart,"-",iend)
    println(f,"#SBATCH --mem=8gb")
    println(f,"#SBATCH --time=",walltime)
    println(f,"#SBATCH --output=\"/scratch2/vdwivedi/disorder/julia/runs/$projname/logs/run%a.log\"")
    println(f,"#SBATCH --job-name=\"$projname\"")
    println(f,"#SBATCH --mail-user=vdwivedi@uni-koeln.de")
    println(f,"#SBATCH --mail-type=BEGIN,END\n")
    println(f,"echo \"Starting. The current time is \" \$(date)")
    println(f,"module load julia")
    println(f,"echo -e \"Loaded Julia module\\n\"\n")
    println(f,"cd /scratch2/vdwivedi/disorder/julia/runs/$projname")
    println(f,"julia  $codepath \${SLURM_ARRAY_TASK_ID}")
    println(f,"echo -e \"\\nDone. The current time is \" \$(date)")
end


# =================================================================

# Create the project directory and the file for the list of parameter values
proj_dir = string(pwd(),"/",projname)
mkpath(string(proj_dir,"/logs"))
println("Created the project directory  ", proj_dir)

jobfile = string(proj_dir,"/job_list.txt")
println("Writing the jobs at ", jobfile)
fjobs = open(jobfile, "w")
println(fjobs,"jobID\tJ0\tJdelta\tJz\tW\tE\tpbc\tLy\tscale")

# Iterator to take all combinations of variables
using Base.Iterators
iter = product(J0list, Jdeltalist, Jzlist, Wlist, Elist, pbclist)
jobsize = length(Wlist)*length(Elist)*length(J0list)*length(Jdeltalist)*length(Jzlist)*length(pbclist)
ord = 1:6 # Reordering vector for the iterator to write the parameters in the correct order

# Need this let environment to be able to update iscr and jobID within the nested loops
let
	jobID = 0
	iscr = 0
	for Ly in Lylist
		for scale in sclist
			# Write out the Bash script
			if bashflag==true
				walltime = computetime(Ly, scale)
				istart = iscr*jobsize + 1
				iend = (iscr+1)*jobsize
				iscr += 1

				# Generate the script file 
				scrname = string(proj_dir,"/jscr",iscr,".sh")
				println("Generating ", scrname)
				scrfile = open(scrname, "w")
				genscript(scrfile, istart, iend, walltime)
				close(scrfile)
			end

			# Write out the parameter values to the joblist file 
			for x in iter
				jobID += 1
				str = ""     # String to store a single row 
				for i in ord
					str *= "$(x[i])\t"
				end
			    println(fjobs,"$jobID\t$str$Ly\t$scale")
			end
		end
	end
end


close(fjobs)
println("Done")
