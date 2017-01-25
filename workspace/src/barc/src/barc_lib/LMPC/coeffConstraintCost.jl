
function coeffConstraintCost(mpcTraj::MpcTrajectory, mpcCoeff::MpcCoeff, posInfo::PosInfo, mpcParams::MpcParams,lapStatus::LapStatus)
    # this computes the coefficients for the cost and constraints

    # Outputs: 
    # coeffConst
    # coeffCost

    # Read Inputs
    s               = posInfo.s
    s_target        = posInfo.s_target

    # Parameters
    Order           = mpcCoeff.order                # interpolation order for cost and constraints
    pLength         = mpcCoeff.pLength              # interpolation length for polynomials
    delay_df        = mpcParams.delay_df

    selected_laps = zeros(Int64,2)
    selected_laps[1] = lapStatus.currentLap-1                                   # use previous lap
    selected_laps[2] = lapStatus.currentLap-2                                   # and the one before
    if lapStatus.currentLap >= 4
        selected_laps[2] = indmin(mpcTraj.cost[1,2:lapStatus.currentLap-2])+1      # and the best from all previous laps
    end

    # Select the old data
    oldS            = mpcTraj.closedLoopSEY[:,1,selected_laps]::Array{Float64,3}
    oldeY           = mpcTraj.closedLoopSEY[:,2,selected_laps]::Array{Float64,3}
    oldePsi         = mpcTraj.closedLoopSEY[:,3,selected_laps]::Array{Float64,3}
    oldV            = mpcTraj.closedLoopSEY[:,4,selected_laps]::Array{Float64,3}
    oldRho          = mpcTraj.closedLoopSEY[:,5,selected_laps]::Array{Float64,3} 

    

    # TODO: Check if using s instead of s_total causes problems
    # Compute the index
    idx_s = zeros(Int64,2)
    for i = 1:2
        DistS = ( s - oldS[:,i] ).^2
        idx_s[i] = Int64(findmin(DistS)[2])
    end              
  

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
            costVector = mpcTraj.cost[vec_range[i],i]                               # decreases in equal steps
            mpcCoeff.coeffCost[:,i] = MatrixInterp[:,:,i]\costVector           # interpolate this vector with the given s
    end
    nothing
end

