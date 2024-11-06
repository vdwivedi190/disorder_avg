# This script generates the list of parameter values for which the disorder-averaging is to be performed and 
# puts it in a directory created for the project, along with subdirectories for run logs. 
# It also optionally generates the bash scripts (one for each Ly) required to run it at a computing cluster. It was written for the CHEOPS cluster at the University of Cologne. 

# PARAMETERS
# ============
# The project name is also the name of the subdirectory that will be created. 
# The relevant system parameters here are the mass (m), Fermi energy (E), and the disorder strength (W)
# One can also scan over a range of the width (Ly) and the aspect ratio (scl = Lx/Ly)
# The type of disoder can be chosen by setting dtype (described below)
# The variable pbc (= 0 or 1) sets the periodic boundary condition

projname = "test"

Lylist = [10 14 20 28 40 56]     # Approximate powers of sqrt(2) times 10
sclist = [1]                # scales = Lx/Ly

mlist = [-0.1 0.1]
Wlist = 0.0:1.0:4.0        # Disorder strength
Elist = [0.0]              # Energy   
pbclist = [1]


# Type of disorder:
#    1,2 = diagonal,sigma_z disorder on every site
#    3,4 = diagonal,sigma_z disorder with density rho, i.e, on rho*Ly unit cells chosen randomly for each iteration

dtype = 2
rholist = [0.1]             # Disorder density (only relevant for dtype = 3,4)


bashflag = false  # Generate bash scripts for HPC cluster if true 

# =================================================================
# The following are relevant only for the bash script
ncpus = 4   
codepath = "./runCI_cluster.jl"

# Estimate the time required and format to a string suitable for the CHEOPS (Cologne HPC cluster) bash script
function computetime(Ly, scale)
    runtime = 2e-6 * Ly^3 * scale + 900 # In seconds
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


# Function to write out the bash scripts for CHEOPS
function genscript(f,istart,iend,walltime)
    println(f,"#!/bin/bash -l\n")
    println(f,"#SBATCH --ntasks=1")
    println(f,"#SBATCH --cpus-per-task=$ncpus")
    println(f,"#SBATCH --array=",istart,"-",iend)
    println(f,"#SBATCH --mem=1gb")
    println(f,"#SBATCH --time=",walltime)
    println(f,"#SBATCH --output=\"/scratch2/vdwivedi/disorder/julia/runs/$projname/logs/run%a.log\"")
    println(f,"#SBATCH --job-name=\"$projname\"")
    println(f,"#SBATCH --mail-user=vdwivedi@uni-koeln.de")
    println(f,"#SBATCH --mail-type=BEGIN,END\n")
    println(f,"echo \"Starting. The current time is \" \$(date)")
    println(f,"module load julia")
    println(f,"echo -e \"Loaded Julia module\\n\"\n")
    println(f,"cd /scratch2/vdwivedi/disorder/julia/runs/$projname")
    println(f,"julia  $codepath \${SLURM_ARRAY_TASK_ID} $dtype")
    println(f,"echo -e \"\\nDone. The current time is \" \$(date)")
end


# =================================================================

# Create the project directory and the file for the list of parameter values
proj_dir = string(pwd(),"/",projname)
mkpath(string(proj_dir,"/logs"))
println("Created the project directory  ", proj_dir)

jobfile = string(proj_dir,"/job_list.txt")
println("Writing the jobs at ", jobfile)
fjobs = open(string(proj_dir,"/job_list.txt"), "w")


# Iterator to take all combinations of variables
using Base.Iterators
if dtype >= 3
	iter = product(Wlist, rholist, Elist, mlist, pbclist)
	#  Note that the definition of jobsize does not include length(Lylist) or length(scllist)
	jobsize = length(Wlist)*length(rholist)*length(Elist)*length(mlist)*length(pbclist)
	ord = [4 1 2 3 5] # Reordering vector for the iterator to write the parameters in the correct order
	println(fjobs,"jobID\tm\tW\trho\tE\tpbc\tLy\tscale")
else
	iter = product(Wlist, Elist, mlist, pbclist)
	#  Note that the definition of jobsize does not include length(Lylist) or length(scllist)
	jobsize = length(Wlist)*length(Elist)*length(mlist)*length(pbclist)	
	ord = [3 1 2 4] # Reordering vector for the iterator to write the parameters in the correct order
	println(fjobs,"jobID\tm\tW\tE\tpbc\tLy\tscale")
end


# Need this let environment to be able to update iscr and jobID within the nested loops
let
	jobID = 0
	iscr = 0
	for Ly in Lylist
		for scale in sclist
			if bashflag==true
				walltime = computetime(Ly, scale)
				#  Compute start and end indices for the script 
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
				str = ""  # String to store a single row 
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
