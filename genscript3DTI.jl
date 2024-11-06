# This script generates the list of parameter values for which the disorder-averaging is to be performed and 
# puts it in a directory created for the project, along with subdirectories for run logs. 
# It also optionally generates the bash scripts (one for each Ly) required to run it at a computing cluster. It was written for the CHEOPS cluster at the University of Cologne. 

# PARAMETERS
# ============
# The project name is also the name of the subdirectory that will be created. 
# The relevant system parameters here are the mass (m), Fermi energy (E), net flux (Phi) and the disorder strength (W)
# One can also scan over a range of the width (Ly=Lz), width of the boundary (Lb), and the aspect ratio (scl = Lx/Ly)
# The type of disoder can be chosen by setting dtype (described below)
# The variable pbc (= 0 or 1) sets the periodic boundary condition

projname = "test_3dti"

# Each system has a square cross section with sides of of Lylist[i] unit cells
Lylist = [16]      # Lz = Ly
Lblist = [2]       # Width of the boundary
sclist = [1]	   # scale = Lx/(Ly*Lz)

mlist = [0.5]
Wlist = [0.25]	      # Disorder strength
Philist = [0.0, 0.25, 0.5, 1.0]    # Net flux 
Elist = 0.01:0.1:0.51

dtype = 1  # 1 for Anderson disorder at each site, 3 for sparse disorder
rholist = [0.1]  # Disorder density (only relevant if dtype=3)


bashflag = false  # Generate bash scripts for HPC cluster if true 

# =================================================================
# The following are relevant only for the bash script
ncpus = 4   
codepath = "./run3DTI_cluster.jl"

# Estimate the time required and format to a string suitable for the CHEOPS (Cologne HPC cluster) bash script
function computetime(Ly, scale)
    runtime = 3e-7 * (Ly*Ly)^3.85 * scale + 600 # In seconds
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


# ===============================================

# Create the project directory and the file for the list of parameter values
proj_dir = string(pwd(),"/",projname)
mkpath(string(proj_dir,"/logs"))
println("Created the project directory  ", proj_dir)

jobfile = string(proj_dir,"/job_list.txt")
println("Writing the jobs at ", jobfile)
fjobs = open(string(proj_dir,"/job_list.txt"), "w")


proj_dir = string(pwd(),"/",projname)
mkpath(string(proj_dir,"/logs"))
datafile = open(string(proj_dir,"/job_list.txt"), "w")

# Iterator to take all combinations of variables
using Base.Iterators
if dtype >= 3
	iter = product(Wlist, rholist, Elist, mlist, Philist)
	jobsize = length(Wlist)*length(rholist)*length(Elist)*length(mlist)*length(Philist)
	ord = [4 1 2 3 5] # Reordering vector for the iterator to write the parameters in the correct order
	println(fjobs,"jobID\tm\tW\trho\tE\tPhi\tLy\tscale")
else
	iter = product(Elist, Wlist, mlist, Philist)
	jobsize = length(Wlist)*length(Elist)*length(mlist)*length(Philist)
	ord = [3 2 1 4] # Reordering vector for the iterator to write the parameters in the correct order
	println(fjobs,"jobID\tm\tW\tE\tPhi\tLy\tLz\tLb\tscale")
end

jobsize *= length(Lblist)

# Need this let environment to be able to update iscr and jobID within the nested loops
let
	jobID = 0
	iscr = 0
	for Ly in Lylist
		for scale in sclist
			# First write the script
			walltime = computetime(Ly, scale)
			istart = iscr*jobsize + 1
			iend = (iscr+1)*jobsize
			iscr += 1

		    scrname = string(proj_dir,"/jscr",iscr,".sh")
			println("Generating ", scrname)
			scrfile = open(scrname, "w")
		    genscript(scrfile, istart, iend, walltime)
		    close(scrfile)

			# Now write out the data points
			for Lb in Lblist
				for x in iter
					jobID += 1
					str = ""
					for i in ord
						str *= "$(x[i])\t"
					end
					println(fjobs,"$jobID\t$str$Ly\t$Ly\t$Lb\t$scale")
				end
			end
		end
	end
end


close(fjobs)
println("Done")
