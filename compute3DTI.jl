@everywhere begin 

# Indices i,j run along axes y,z, respectively.
@everywhere index(i,j) = LinearIndices((Ly,Lz))[i,j]
@everywhere site(i) = CartesianIndices((Ly,Lz))[i]

# Function to construct the J and M matrices for the 3dTI nanowires (without magnetic field)
@everywhere function construct_mat(m,Ly,Lz,pbc)

    # Pauli matrices
    sig_0 =[1 0; 0 1]
    sig_x= [0 1;1 0]
    sig_y= [0 -im; im 0 ]
    sig_z =[1 0; 0 -1]

    # Tensor products to get a set of four anticommuting 4x4 matrices
    alpha_x = stdform(kron(sig_z,sig_x))
    alpha_y = stdform(kron(sig_z,sig_y))
    alpha_z = stdform(kron(sig_z,sig_z))
    beta = stdform(kron(sig_x,sig_0))

    # Coefficients of 1, e^{i k_{x,y,z}} in H(k)
    Mmat = (3-m) * beta
    Jxmat = (im*alpha_x + beta)/2
    Jymat = (im*alpha_y + beta)/2
    Jzmat = (im*alpha_z + beta)/2

    Diag = Diagonal(ones(Ly*Lz))
    J = kron(Diag,Jxmat)  # The transfer marix works along the x-axis
    M = kron(Diag,Mmat)   # Purely on-site terms

    for i in 1:(Ly-1)
        for j in 1:(Lz-1)
            # n1 and n2 represents the indices of two nearest-neighbor sites
            n1 = index(i,j)

            # Bonds along the y-axis 
            n2 = index(i+1,j)
            M[uc*(n1-1)+1:uc*n1, uc*(n2-1)+1:uc*n2] .= Jymat
            M[uc*(n2-1)+1:uc*n2, uc*(n1-1)+1:uc*n1] .= Jymat'

            # Bonds along the z-axis 
            n2 = index(i,j+1)
            M[uc*(n1-1)+1:uc*n1, uc*(n2-1)+1:uc*n2] .= Jzmat
            M[uc*(n2-1)+1:uc*n2, uc*(n1-1)+1:uc*n1] .= Jzmat'
        end
    end

    return(J,M)
end

# Function to construct the J and M matrices for the 3dTI in presence of a magnetic field
# (Using the model from Cook and Franz, PRB 84, 201105(R), (2011): Eq (9) 
@everywhere function construct_mat_flux(m,Ly,Lz,Phi)

    # Pauli matrices
    sig_x= [0 1;1 0]
    sig_y= [0 -im; im 0 ]
    sig_z =[1 0; 0 -1]
    sig_0 =[1 0; 0 1]

    # Tensor to get a set of anticommuting 4x4 matrices
    alpha_x = stdform(kron(sig_y,sig_0))
    alpha_y = stdform(kron(sig_z,sig_z))
    alpha_z = -stdform(kron(sig_z,sig_y))
    beta = stdform(kron(sig_x,sig_0))

    # Coefficients of 1, e^{i k_{x,y,z}} in H(k)
    Mmat = (3-m) * beta
    Jxmat = (im*alpha_x + beta)/2
    Jymat = (im*alpha_y + beta)/2
    Jzmat = (im*alpha_z + beta)/2
    
    Diag = Diagonal(ones(Ly*Lz))
    J = kron(Diag,Jxmat)
    M = kron(Diag,Mmat)

    B = Phi/((Ly-1)*(Lz-1))  # Magnetic field

    for i in 1:Ly
        for j in 1:Lz
            # n1 and n2 represents the indices of two nearest-neighbor sites
            n1 = index(i,j)
            
            # Bonds along the y-axis 
            if i < Ly
                n2 = index(i+1,j)
                M[uc*(n1-1)+1:uc*n1, uc*(n2-1)+1:uc*n2] .= Jymat
                M[uc*(n2-1)+1:uc*n2, uc*(n1-1)+1:uc*n1] .= Jymat'
            end
        
            # Bonds along the z-axis (Landau gauge with Az(y) ≠ 0 )
            if j < Lz
                n2 = index(i,j+1)
                M[uc*(n1-1)+1:uc*n1, uc*(n2-1)+1:uc*n2] .= Jzmat * exp(2*pi*im*i*B)
                M[uc*(n2-1)+1:uc*n2, uc*(n1-1)+1:uc*n1] .= Jzmat' * exp(-2*pi*im*i*B)
            end

        end
    end
    return(J,M)
end

# A straight line profile for a smooth boundary of thickness xi 
@everywhere function wire_profile(y::Float64, xi::Float64)
    if y <= xi
        return 0
    elseif y >= 1-xi
        return 1
    else
        return (y-xi)/(1-2*xi)
    end
end



# This function Generate an array with elements 0 or 1 denoting whether a given site has disorder or not. 
# The resulting on-site matrix contains Wd[2] disordered sites on average.
@everywhere function add_disorder!(disordered_M::Array{Complex{Float64},2},Wd::Array{Float64,1},dmat::Array{Complex{Float64},2},Nuc::Int64,iter::Int64)    
    
    if Wd[1] > 0   # Only modify disordered_M if the disorder strength is nonzero.
        if Wd[2] == 1.0    
            sites = ones(Nuc)   # Disorder density = 1, i.e., each site has a random on-site chemical potential
        else
            # Disorder density < 1. 
            # The array sites contains a random list of sites (with the given disorder density) on which a disorder is added. 
            sites = floor.(rand(Nuc) .+ Wd[2])        
        end

        # List of random on-site potentials for the disordered sites 
        rlist = rand(Nuc)
        rlist .= Wd[1] .* (rlist .- 0.5)

        @inbounds @simd for j = 1:Nuc
            disordered_M[uc*(j-1)+1:uc*j, uc*(j-1)+1:uc*j] .+= sites[j] * rlist[j] * dmat
        end
    end

    # Add the large potential at the boundary to change the cross-section of the wire
    U0 = 20;    # Boundary potential (to choose a smooth boundary)
    xi = 0.0;   # Smoothness of the boundary     
    disordered_M .+= U0 * wire_profile(iter/Lx,xi) * Bmask;
end



@everywhere function perform_job(my_job, dtype, qrflag)

    # Number of elements specifying the disorder: 1 (disorder strength) for disorder on each site and 2 (disorder strength and density) for sparse disorder
    dlen = dtype < 3 ? 1 : 2

    # =================================================
    # Parse the job details
    global jobID = Int64.(my_job[1])
    m = my_job[2]
    W = dlen == 1 ? [my_job[3], 1.0] : [my_job[3], my_job[4]]    
    E = my_job[3+dlen]
    Phi = my_job[4+dlen]
    global Ly = Int64.(my_job[5+dlen])
    global Lz = Int64.(my_job[6+dlen])
    Lb = Int64.(my_job[7+dlen])
    scale = Int64.(my_job[8+dlen])
    global Lx = scale*Ly*Lz

    if 2*Lb > min(Ly,Lz)
        Lb = floor(min(Ly,Lz)/2)
        println("Given boundary width larger than the cross section. Setting it to ", Lb, ".")
    end

    print("Imported data. Computing the boundary sites...")    
    
    # =================================================
	# Define a mask demarcating the boundary
        
    # The matrix boundary is 1 at the "boundary sites" and zero elsewhere
    # This section is relevant only when Lb > 0. 
    boundary = zeros(Ly*Lz)
    for b in 1:Lb
        # Set the b'th element from the boundary along the y-axis
        for i in 1:Ly
            boundary[index(i,b)] = 1.0
            boundary[index(i,Lz-b+1)] = 1.0
        end

        # Set b'th element from the boundary along the z-axis
        for i in 1:Lz
            boundary[index(b,i)] = 1.0
            boundary[index(Ly-b+1,i)] = 1.0
        end
    end
    println("Done")

    println("Mask for the cross section:")
    display(Int64.(boundary[index(:,:)]))
    println()

    # Construct a diagonal matrix with identity matrices at each boundary site (Anderson disorder).
    # Can add a boundary potential to the Hamiltonian by adding U*Bmask
    dmat = stdform(Diagonal(ones(uc)))
    global Bmask = kron(diagm(0 => boundary), dmat)

    # =================================================
    # Construct J, M and call getLyapunov() 
    worker_id = myid()
    host_id = gethostname()
    start_time= now()
    println("Starting my job $(jobID) on $(host_id) at $(start_time) ")
    
    J,M = construct_mat_flux(m,Ly,Lz,Phi)
    q = 3  # Number of steps between QR decompositions
    dmat = stdform(dmat)

    t_iter = 3e-7 * (Ly*Ly)^2.85    # Estimatated emperically from previous runs
    blocksize = Int64(1000*ceil(10/t_iter))  # Dump data every ~3 hours

    # calculate the Lyapunov's for the given job data
    λ_list,Q_prev,R = getLyapunovList(J,M,E,Ly,Lx,W,dmat,q,qrflag,blocksize)

    finish_time= now()
    time_taken = Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(finish_time) - Dates.DateTime(start_time)))
	println("Finished my job $(jobID) on $(host_id) at time $(finish_time) ")

	# Update the log file
    open("./run_log.txt", "a") do f
        write(f,"$(jobID)\t $m\t $(W[1])\t $Lx\t $Ly\t $Lz\t $Lb\t $worker_id\t\t $host_id\t $time_taken\n")
	end
end

end