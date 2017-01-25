    include("../workspace/src/barc/src/barc_lib/classes.jl")
    include("../workspace/src/barc/src/barc_lib/LMPC/functions.jl")
    include("../workspace/src/barc/src/barc_lib/LMPC/coeffConstraintCost.jl")


    mpcCoeff                    = MpcCoeff() #m: coefficients for SS-traj and cost approximation
    mpcTraj                     = MpcTrajectory()   # data type to store the mpc trajectories and cost

    buffersize = 10000
    mpcTraj.closedLoopSEY       = zeros(buffersize,7,30)
    mpcTraj.inputHistory        = zeros(buffersize,3,30)
    mpcTraj.cost                = zeros(buffersize,30)
    mpcTraj.idx_end             = zeros(30) # all data points until between 0 <= s <= lapLength
    mpcTraj.count               = ones(30) #total number of saved data points for each lap

    mpcTraj.closedLoopSEY[1:buffersize,1,1] = linspace(0,19.11,buffersize)
    mpcTraj.closedLoopSEY[1:buffersize,1,2] = linspace(0,19.11,buffersize)

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(mpcCoeff.order+1,2)
    mpcCoeff.coeffConst         = zeros(mpcCoeff.order+1,2,4)
    mpcCoeff.pLength            = 2*12  

    s               = 19.11/2.0
    s_target        = 19.11

    # Parameters
    Order           = mpcCoeff.order                # interpolation order for cost and constraints
    pLength         = mpcCoeff.pLength              # interpolation length for polynomials

    selected_laps = [1,2]

    # Select the old data
    oldS            = mpcTraj.closedLoopSEY[:,1,selected_laps]::Array{Float64,3}
    oldeY           = mpcTraj.closedLoopSEY[:,2,selected_laps]::Array{Float64,3}
    oldePsi         = mpcTraj.closedLoopSEY[:,3,selected_laps]::Array{Float64,3}
    oldV            = mpcTraj.closedLoopSEY[:,4,selected_laps]::Array{Float64,3}
    oldRho          = mpcTraj.closedLoopSEY[:,5,selected_laps]::Array{Float64,3} 

    

    # TODO: Check if using s instead of s_total causes problems
    # Compute the index
    DistS = ( s - oldS[:,1] ).^2
    DistS2 = ( s - oldS[:,2] ).^2

    idx_s = zeros(Int64,2)
    idx_s[1] = Int64(findmin(DistS)[2])              
    idx_s[2] = Int64(findmin(DistS2)[2])          

    vec_range = (idx_s[1]:idx_s[1]+pLength,idx_s[2]:idx_s[2]+pLength)
    # Create the vectors used for the interpolation
    # bS_vector contains the s-values for later interpolation
    bS_Vector       = zeros(pLength+1,2)
    for i=1:pLength+1
        bS_Vector[i,1] = oldS[vec_range[1][i]]
        bS_Vector[i,2] = oldS[vec_range[2][i]]
    end


    # Create the Matrices for the interpolation
    MatrixInterp = zeros(pLength+1,Order+1,2)

    for k = 0:Order
        MatrixInterp[:,Order+1-k,:] = bS_Vector[:,:].^k
    end
    
    # Compute the coefficients
    for i=1:2
        mpcCoeff.coeffConst[:,i,1]    = MatrixInterp[:,:,i]\oldeY[vec_range[i]]
        mpcCoeff.coeffConst[:,i,2]    = MatrixInterp[:,:,i]\oldePsi[vec_range[i]]
        mpcCoeff.coeffConst[:,i,3]    = MatrixInterp[:,:,i]\oldV[vec_range[i]]
        mpcCoeff.coeffConst[:,i,4]    = MatrixInterp[:,:,i]\oldRho[vec_range[i]]
    end

    # Finished with calculating the constraint coefficients
    
    # Now compute the final cost coefficients
    # The Q-function contains for every point in the sampled safe set the minimum cost-to-go-value
    # These values are calculated for both old trajectories

    for i=1:2 
            @show(vec_range[i])  
            costVector = mpcTraj.cost[vec_range[i],i]                               # decreases in equal steps
            @show(costVector)
            mpcCoeff.coeffCost[:,i] = MatrixInterp[:,:,i]\costVector           # interpolate this vector with the given s
    end