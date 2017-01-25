
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

    # Redo component calculation for the two selected_laps
    for kkk = 1:2
        lapNum = selected_laps[kkk]
        oldS            = mpcTraj.closedLoopSEY[:,1,lapNum]
        oldeY           = mpcTraj.closedLoopSEY[:,2,lapNum]
        oldePsi         = mpcTraj.closedLoopSEY[:,3,lapNum]
        oldV            = mpcTraj.closedLoopSEY[:,4,lapNum]
        oldRho          = mpcTraj.closedLoopSEY[:,5,lapNum]

    

        # TODO: Check if using s instead of s_total causes problems
        # Compute the index
        DistS = ( s - oldS ).^2
        idx_s = Int64(findmin(DistS)[2])
        
        vec_range = idx_s:idx_s+pLength
        bS_Vector = oldS[vec_range]


        # Create the Matrices for the interpolation
        MatrixInterp = zeros(pLength+1,Order+1)

        for k = 0:Order
            MatrixInterp[:,Order+1-k] = bS_Vector.^k
        end
        
        # Compute the coefficients
        mpcCoeff.coeffConst[:,kkk,1]    = MatrixInterp[:,:]\oldeY[vec_range]
        mpcCoeff.coeffConst[:,kkk,2]    = MatrixInterp[:,:]\oldePsi[vec_range]
        mpcCoeff.coeffConst[:,kkk,3]    = MatrixInterp[:,:]\oldV[vec_range]
        mpcCoeff.coeffConst[:,kkk,4]    = MatrixInterp[:,:]\oldRho[vec_range]

        # Finished with calculating the constraint coefficients
        
        # Now compute the final cost coefficients
        # The Q-function contains for every point in the sampled safe set the minimum cost-to-go-value
        # These values are calculated for both old trajectories

        costVector = mpcTraj.cost[vec_range,lapNum]                               # decreases in equal steps
        mpcCoeff.coeffCost[:,kkk] = MatrixInterp[:,:]\costVector           # interpolate this vector with the given s
    end
    nothing
end

