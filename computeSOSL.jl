@everywhere begin 
    
# Function to construct the J and M matrices for the SOSL model on the Shastry-Sutherland lattice
# See Physical Review Research 2 (4), 043159 (2020) for details on the model
@everywhere function construct_mat(J0,Jdelta,Jdiag,Ly,pbc)
	J1 = J0+Jdelta
	J2 = J0-Jdelta

    Diag = Diagonal(ones(Ly))
    Udiag = diagm(1 => ones(Ly-1))
    Ldiag = diagm(-1 => ones(Ly-1))

    # Define the hopping matrices for each 4-site unit cell along x and y 
    Jxmat1 = im*J2*[0 0 0 0; 1 0 0 0; 0 0 0 1; 0 0 0 0]
    Jxmat2 = im*Jdiag*[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 0 0]
	Jymat = im*J2*[0 0 0 0; 0 0 0 0; 0 1 0 0; -1 0 0 0]

    # Define the intra-unit-cell bonds along sides and diagonals of the unit cell
    Mmat1 = im*J1*[0 1 0 -1; -1 0 1 0; 0 -1 0 -1; 1 0 1 0]
    Mmat2 = im*Jdiag*[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 -1 0 0]

    # Construct J and M for a strip of width Ly unit cells (so that J,M are 4Ly x 4Ly matrices)
    J = kron(Diag,Jxmat1) + kron(Udiag,Jxmat2)
	M = kron(Diag,Mmat1+Mmat2) + kron(Udiag,Jymat) + kron(Ldiag,adjoint(Jymat))

    if pbc == 1
        LdiagPBC = diagm((Ly-1) => ones(1))
        UdiagPBC = diagm(-(Ly-1) => ones(1))
        M += + kron(UdiagPBC,Jymat) + kron(LdiagPBC,adjoint(Jymat))
        J += kron(UdiagPBC,Jxmat2)
    end

    # Real space Hamiltonian: Can write to a file for exact diagonalization
    # H = kron(Diag,M) + kron(Udiag,J) + kron(Ldiag,adjoint(J))
    # writedlm("./H.dat",H)

    return(J,M)
end

# Bonds on which the sign of the hopping must be flipped to generate four possible vison configurations
# The entries are site indices (1-4), and 6 represents site #2 in the next unit cell, which site #3 in a unit cell connects to. 
@everywhere bonds_to_flip = [1 2; 2 3; 3 4; 3 6]


# Function to add the flux-disorder to the 4-plaquettes
# The disorder strength Wd is 1-element array, whose only element is the probability of a vison
@everywhere function add_disorder!(disordered_M::Array{Complex{Float64},2},Wd::Array{Float64,1},dmat::Array{Complex{Float64},2},Nuc::Int64,iter::Int64)
	# Generate a random matrix with entries ±1 with probability 1-w and w
    rmat = sign.(rand(Nuc,uc) .- Wd[1])

	N = Nuc*uc

	for j = 1:Nuc              # Loop over the unit cells along y
		for k = 1:uc           # Loop over the sites within a unit cell 
			if rmat[j,k] == -1
				# Compute indices corresponding to the two sites connected by the bond 
				c1 = uc*(j-1) + bonds_to_flip[k,1]
				c2 = uc*(j-1) + bonds_to_flip[k,2]

                # Wrap around (assumes periodic BC - trivial gauge transformation for open BC)
				c2 <= N ? c2 : c2 -= N     

				# Flip the sign of the requisite bond
				disordered_M[c1,c2] *= -1
				disordered_M[c2,c1] *= -1
			end
		end
	end
end

# Function to parse a line of parameter values and invoke getLyapunov() 
@everywhere function perform_job(my_job,qrflag)

	# Parse the job details
    global jobID = Int64.(my_job[1])
    J0 = my_job[2]
    Jdelta = my_job[3]
    Jdiag = my_job[4]
    W = [my_job[5]]   # The input W is expected to be a vector
    E = my_job[6]
    pbc = my_job[7]
    Ly = Int64.(my_job[8])
    scale = Int64.(my_job[9])
    Nx = scale*Ly

	worker_id = myid()
    host_id = gethostname()
    start_time= now()
    println("Starting my job $(jobID) on $(host_id) at time $(start_time) ")

	J,M = construct_mat(J0,Jdelta,Jdiag,Ly,pbc)
	q = 1 # Number of steps between QR decompositions

	# The following matrix is needed only for diagonal disorders, so its value is irrelevant for the present computation
	dmat = stdform([1 0; 0 1])

	t_iter = 1e-6 * Ly^2.6
	blocksize = Int64(1000*ceil(10/t_iter))  # Dump data every ~3 hours

	# Calculate the Lyapunov's for the given job data
	λ_list,Q_prev,R = getLyapunovList(J,M,E,Ly,Nx,W,dmat,q,qrflag,blocksize)

    finish_time = now()
    time_taken = Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(finish_time) - Dates.DateTime(start_time)))
    println("Finished my job $(jobID) on $(host_id) at time $(finish_time) ")    
    
    # Update the log file
	open("./run_log.txt", "a") do f
         write(f,"$(jobID)\t$(J0)\t$(Jdelta)\t$(Jdiag)\t$W\t$E\t$(Ly)\t$(Nx)\t$(worker_id)\t$host_id\t $time_taken\n")
    end

end

end