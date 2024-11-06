@everywhere begin 
	
# Construct the J and M matrices for Chern insulator 
# (See https://doi.org/10.1103/PhysRevB.93.134304 for their definitions)
@everywhere function construct_mat(m,Ly,pbc)
	Diag=Diagonal(ones(Ly))
	Udiag = diagm(1 => ones(Ly-1))
    Ldiag = diagm(-1 => ones(Ly-1))

	sig_x= [0 1;1 0]
	sig_y= [0 -im; im 0 ]
	sig_z =[1 0; 0 -1]

	Jxmat = -(im/2 ) * (sig_x - 1im*sig_z)
	Jymat = im/2 * (sig_y + 1im* sig_z)
	Mmat = (2-m) * sig_z

    J = kron(Diag,Jxmat)
	M = kron(Diag,Mmat) + kron(Udiag,Jymat) + kron(Ldiag,Jymat')
	
	if pbc == 1
        LdiagPBC = diagm((Ly-1) => ones(1))
        UdiagPBC = diagm(-(Ly-1) => ones(1))
        M += + kron(UdiagPBC,Jymat) + kron(LdiagPBC,Jymat')
        J += kron(UdiagPBC,Jxmat)
    end

    return(J,M)
end


# Generate the on-site disorder for the Chern insulator
@everywhere function add_disorder!(disordered_M::Array{Complex{Float64},2},Wd::Array{Float64,1},dmat::Array{Complex{Float64},2},Nuc,iter)
    # Generate an array with elements 0 or 1 denoting whether a given site has disorder or not. The array contains Wd[2] disordered sites on average

	if Wd[1] > 0   # Do not modify disordered_M if the disorder strength is zero.
    	if Wd[2] == 1.0
        	sites = ones(Nuc)
	    else
    	    sites = floor.(rand(Nuc) .+ Wd[2])        # Disorder with density < 1
	    end

	    # Generate list of random potentials
    	rlist = rand(Nuc)
	    rlist .= Wd[1] .* (rlist .- 0.5)

		@inbounds @simd for j = 1:Nuc
    	    disordered_M[uc*(j-1)+1:uc*j, uc*(j-1)+1:uc*j] .+= sites[j] * rlist[j] * dmat
	    end

		# The singularity of a single block sometimes makes M singular
		# Add a small off diagonal bit to mitigate that issue
		rlist2 = rand(Nuc*uc)
        rlist2 .= 1e-6 * (rlist2 .- 0.5)
        disordered_M .+= Diagonal(rlist2)
	end
end

# Function to parse a line of parameter values and invoke getLyapunov() 
@everywhere function perform_job(my_job, dtype, qrflag)	

	# Number of elements specifying the disorder: 1 (disorder strength) for disorder on each site and 2 (disorder strength and density) for sparse disorder
	dlen = dtype < 3 ? 1 : 2

	# Parse the job details
    global jobID = Int64.(my_job[1])
	m = my_job[2]
	W = (dlen == 1) ? [my_job[3], 1.0] : [my_job[3], my_job[4]]
	E = my_job[3+dlen]
    pbc = my_job[4+dlen]
	Ly = Int64.(my_job[5+dlen])
	scale = Int64.(my_job[6+dlen])
	Nx = scale*Ly

    # =================================================
	# Define the disorder matrix based on the disorder type
	
	if dtype == 1
		println("Disorder type: Diagonal disorder at each site")
		dmat=[1 0; 0 1]
	elseif dtype == 2
		println("Disorder type: Sigma_z disorder at each site")
		dmat=[1 0; 0 -1]
	elseif dtype == 3
		println("Disorder type: Diagonal disorder at density ", W[2])
		dmat=[1 0; 0 1]
	elseif dtype == 4
		println("Disorder type: Sigma_z disorder at density ", W[2])
		dmat=[1 0; 0 -1]
	else
		println("Invalid disorder specification. Aborting!")
		exit()
	end	

    # =================================================
    # Construct J, M and call getLyapunov() 
	worker_id = myid()
    host_id = gethostname()
	start_time= now()
	println("Starting my job number $(jobID) on $(host_id) at $(start_time) ")
	println("The current directory is $(pwd()) ")
		
	J,M = construct_mat(m,Ly,pbc)
	q = 3  # Number of steps between QR decompositions
	dmat = stdform(dmat)

	t_iter = 3e-7 * Ly^2.6   # Estimatated emperically from previous runs
    blocksize = Int64(1000*ceil(10/t_iter))  # Dump data every ~3 hours
		
	# calculate the Lyapunov's for the given job data
	λ_list,Q_prev,R = getLyapunovList(J,M,E,Ly,Nx,W,dmat,q,qrflag,blocksize)

    finish_time = now()
    time_taken = Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(finish_time) - Dates.DateTime(start_time)))
	println("Finished my job $(jobID) on $(host_id) at time $(finish_time) ")

    
	# Update the log file
    open("./run_log.txt", "a") do f
		write(f,"$(jobID)\t$m\t$(W[1])\t$(Ly)\t $(Nx)\t$(worker_id)\t$host_id\t\t$time_taken\n")
	end

end


end