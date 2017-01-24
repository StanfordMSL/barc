
function coeffConstraintCost(oldTraj::OldTrajectory, mpcCoeff::MpcCoeff, posInfo::PosInfo, mpcParams::MpcParams,lapStatus::LapStatus)
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
    if lapStatus.currentLap >= 5
        selected_laps[2] = indmin(mpcTraj.cost[1,2:lapStatus.currentLap-2])+1      # and the best from all previous laps
    end

    # Select the old data
    oldS            = mpcTraj.closedLoopSEY[:,1,selected_laps]::Array{Float64,3}
    oldeY           = mpcTraj.closedLoopSEY[:,2,selected_laps]::Array{Float64,3}
    oldePsi         = mpcTraj.closedLoopSEY[:,3,selected_laps]::Array{Float64,3}
    oldV            = mpcTraj.closedLoopSEY[:,4,selected_laps]::Array{Float64,3}
    oldRho          = mpcTraj.closedLoopSEY[:,5,selected_laps]::Array{Float64,3} 

    DistS::Array{Float64}   # initialize
    idx_s::Array{Int64}     # initialize
    vec_range::Tuple{UnitRange{Int64},UnitRange{Int64}}
    costVector::Array{Float64}


    # Compute the index
    DistS = ( s - oldS ).^2

    idx_s = findmin(DistS,1)[2]              # contains both indices for the closest distances for both oldS !!

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
            costVector = mpcTraj.cost[vec_range[i]]                               # decreases in equal steps
            mpcCoeff.coeffCost[:,i] = MatrixInterp[:,:,i]\costVector           # interpolate this vector with the given s
    end
    nothing
end

# Notes about oldTrajectory:
# oldTrajectory[1:prebuf] = states before finish line (s < 0)
# oldTrajectory[prebuf+1:prebuf+cost] = states between start and finish line (0 <= s < s_target)
# oldTrajectory[prebuf+cost+1:prebuf+cost+postbuf] = states after finish line (s_target <= s)
# once one lap is over, the states are saved (interval prebuf:s_target)
# during the next lap, the states behind the finish line are appended (interval s_target:postbuf)
