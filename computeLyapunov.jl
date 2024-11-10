# This file contains the function get_LyapunovList() which takes a 2d tight-binding Hamiltonian (described by the hopping matrix and the on-site matrix) and computes the disorder-averaged Lyapunov exponents. 

@everywhere begin

# Standard typecasting for matrices (everything should be a complex matrix)
@everywhere function stdform(mat)
    convert(Array{Complex{Float64},2},Matrix(mat))
end



# In-place submatrix assignments
# Assigns the matrix elements of A to those of B (with dim B < dim A) starting at indices (m,n).
@everywhere function assign_submat!(A::Array{Complex{Float64},2},B::Array{Complex{Float64},2},m,n)
    Blen = size(B,1)
    @inbounds @simd for j in n:(n+Blen-1)
        @inbounds for i in m:(m+Blen-1)
            A[i,j] = B[i-m+1, j-n+1];
        end
    end
end


# Function to compute the disorder-averaged Lyapunov exponents
# ============================================================
# J and M are the hopping and on-site matrices, respectively. 
# E is the energy at which the Lyapunov exponents are computed 
# Ly is the width of the system (number of unit cells)
# Nx is the number of slices over which the disorder-average is performed 
# Wd and dmat denote the strength and matrix-form of the disorder 
# q is the number of steps after a QR-decomposition is performed 
# If qrflag=True, the Q and R matrices are also stored in file. 
# blocksize denotes the number of iterations after which the Lyapunov exponents are saved to a file (to salvage the data in case the program crashes/is terminated before Nx iterations).  

@everywhere function getLyapunovList(J::Array{Complex{Float64},2},M::Array{Complex{Float64},2},E::Float64,Ly,Nx,Wd::Array{Float64,1},dmat::Array{Complex{Float64},2},q,qrflag,blocksize)


#   INITIALIZATIONS
#  ==========================================

    # ===== Parameters (scalars) =========

    # System size
    N = size(M,1)  # Number of degrees of freedom per unit cell
    r = rank(J)
    sizeT = 2*r    # The rank of J determines the size of the transfer matrix

    Nuc = Int(floor(N/uc))  # Number of unit cells

    # Iteration parameters
    Nblocks = Int(floor(Nx/blocksize))
    lastblocksize = Nx - blocksize*Nblocks

    unit = 1.0+0.0*im;

    # Matrix to temporarily store a disordered configuration
    disordered_M = copy(M)

    # Cumulatively store the Lyapunov exponents in λ_list
    λ_list = stdform(zeros(sizeT,1))


    # ===== Constant matrices ========

    # ==== N x N ======
    Id_N = stdform(Diagonal(ones(N)))   # Identity matrix of size: N x N

    # Initialize local Green's function and on-site matrix as the identity matrix 
    G, Mtemp = copy(Id_N), copy(Id_N)


    # ===== r x r =======
    O = stdform(Diagonal(zeros(r)))  # Zero matrix of size: r x r
    Id_r = stdform(Diagonal(ones(r)))  #Identity matrix of size: r X r
    negId = (-1).*Id_r;

    # Subblocks of the on-site Green's function required to compute the transfer matrix
    Gvv, Gvw, Gwv, Gww = copy(Id_r),copy(Id_r),copy(Id_r),copy(Id_r)

    # Copies of above for in-place matrix operations
    Gvv1, Gvw1, Gwv1, Gww1 = copy(Id_r),copy(Id_r),copy(Id_r),copy(Id_r)


    # ===== 2r x 2r =======
    Id_2r = stdform(Diagonal(ones(sizeT)))  #Identity matrix of size 2r X 2r
    T = copy(Id_2r)                         # Initialize the transfer matrix
    A, B = copy(Id_2r), copy(Id_2r)         # Auxilliaries reqd to compute T
    Tcur = copy(Id_2r)                      # Transfer matrix at the current step
    temp = copy(Id_2r)                      # Temporary variable to store products
    Q_prev, R = copy(Id_2r), copy(Id_2r)    # Variables to store QR decompositions


#   LOCAL FUNCTIONS
#  ==========================================

    # Define this to be "global" within a function call to calc_LyapunovList()
    temp_rN = stdform(zeros(r,N))
    
    # In-place multiplier for 3 nonsquare matrices, as needed to compute G_ab
    # Multiples the three (nonsquare) matrices in "list" and overwrites the result on G_ab.
    function comp_overlap!(G_ab, list)
      mul!(temp_rN, list[1], list[2])
      mul!(G_ab, temp_rN, list[3])
    end


    # This performs a single iteration, and stores everything in the variables local to getLyapunovList()
    function iterate!(ind)
        for i = 1:q
            for l = 1:100   # Try up to 100 disorder realizations to get a nonsingular G
                disordered_M .= M;

                # Modified this call on 05.02.20: Added the last argument to this function to
                # include the iteration index (to encode the wire profile for 3dTI wire.)
                add_disorder!(disordered_M,Wd,dmat,Nuc,q*ind+i-1)
                Mtemp .= E.*Id_N .- disordered_M;  

        		F = lu!(Mtemp, check=false)

        		if issuccess(F)
        		    G .= Id_N
        		    ldiv!(F,G)
        		    break
        		end
                println("Encountered a singular disorder configuration at iteration #",ind,"! (l = $l)")
                println("\nTrying a new disorder realization now.\n")
                flush(stdout)
            end

           # Calculate the submatrices needed to assemble the transfer matrix
            comp_overlap!(Gvv,(Vd,G,V))
            comp_overlap!(Gvw,(Wt,G,V))
            comp_overlap!(Gwv,(Vd,G,Wtd))
            comp_overlap!(Gww,(Wt,G,Wtd))

            mul!(Gvv1, Gvv, Xi)
            mul!(Gvw1, Gvw, Xi)
            mul!(Gwv1, Gwv, Xi)
            mul!(Gww1, Gww, Xi)

            assign_submat!(A, Gvv1, 1, 1)
            assign_submat!(A, negId, 1, r+1)
            assign_submat!(A, Gvw1, r+1, 1)
            assign_submat!(A, O, r+1, r+1)

            assign_submat!(B, O, 1, 1)
            assign_submat!(B, Gwv1, 1, r+1)
            assign_submat!(B, negId, r+1, 1)
            assign_submat!(B, Gww1, r+1, r+1)

            # Calculate the transfer matrix using the SVD basis of J
            ldiv!(lu!(A),B);  # Computes A^{-1} B and stores it in B.
            Tcur .= (-1).*B;

            if i==1         # Starting a new q block
                T .= Tcur
            else            # Already inside a q block, take product of all T matrices
                mul!(temp, T, Tcur)
                T .= temp
            end

            # At the end of the q block (possibly with q = 1)
            if i==q
                R .= T*Q_prev   # T ----> T'= T*Q_(x-1), hold the result in R
                F = qr!(R)
                Q_prev = F.Q    #Q_prev stores Q_x for next iteration x+1
            end

          end
    end



#  DEFINE FILE PATHS
# ==============================================================

    # Create directories for outputs
    dir_name = "."
	mkpath(string(dir_name,"/lyaps"))
    ltempfile = string(dir_name,"/lyaps/l_",jobID,".tmp")
    lfile = string(dir_name,"/lyaps/l_",jobID)

    if qrflag == 1
        mkpath(string(dir_name,"/Q"))
        mkpath(string(dir_name,"/R"))

        Qtempfile = string(dir_name,"/Q/Q_",jobID,".tmp")
        Rtempfile = string(dir_name,"/R/R_",jobID,".tmp")
        Qfile = string(dir_name,"/Q/Q_",jobID)
        Rfile = string(dir_name,"/R/R_",jobID)
    end


#   SVD
#  ==========================================

    # Perform reduced SVD of J only once in the beginning (Will need to change that for bond disorders) 
    # The SVD follows the notation J = V.Xi.W', where W' is the transpose conjugate of W
    F = svd(J)
    V = stdform(F.U[:,1:r])
    Wt = stdform(F.Vt[1:r,:])
    Xi = stdform(diagm(0 => F.S[1:r]))

    # Transpose conjugate("dagger") of V and Wt (W transpose) 
    Vd = V'
    Wtd = Wt'



#  THE MAIN LOOP
# ==============================================================

	println("Rank of the hopping matrix = $r")
    # Throw away the first q x Ntrans iterations (transients)
    Ntrans = 10;
    for x in (1:Ntrans)
       iterate!(0)
    end

    println("Done with preliminary iterations.")
    flush(stdout)

    # First Nblocks set of iterations with data dumps at the end
    if Nblocks>0
    	println("Nx = ", Nx,". Temporary data dump every ", blocksize," iterations.")
    	flush(stdout)

        for x in 1:Nblocks
            for y in 1:Int(ceil(blocksize/q))
                iterate!(blocksize*(x-1)+y)
                λ_list .+= log.(abs.( view(R, diagind(R)) )) ./ Nx
            end

            writedlm(ltempfile,real(λ_list), ", ")
            if qrflag == 1
                writedlm(Qtempfile,Q_prev, ", ")
                writedlm(Rtempfile,R, ", ")
            end
            println("Temporary data dump #", x," at ", now())
            flush(stdout)
        end
    end

    # Remaining iterations
    for y in 1:Int(ceil(lastblocksize/q))
        iterate!(blocksize*Nblocks+y)
        λ_list .+= log.(abs.( view(R, diagind(R)) )) ./ Nx
    end

    # Write out the final data
    writedlm(lfile,real(λ_list), ", ")
    if qrflag == 1
        writedlm(Qfile,Q_prev, ", ")
        writedlm(Rfile,R, ", ")
    end
    println("Wrote out the final Lyapunov exponents to ", lfile)
    flush(stdout)

    # Remove the temporary files if they were generated
    if Nblocks>0
        rm(ltempfile)
        if qrflag == 1
            rm(Qtempfile)
            rm(Rtempfile)
        end
    end

    return(λ_list,Q_prev,R)
end


end