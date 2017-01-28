#!/usr/bin/env julia4

using RobotOS
@rosimport barc.msg: ECU, pos_info
@rosimport data_service.msg: TimeData
@rosimport geometry_msgs.msg: Vector3
rostypegen()
using barc.msg
using data_service.msg
using geometry_msgs.msg
using JuMP
using Ipopt
using JLD

# log msg
include("barc_lib/classes.jl")
include("barc_lib/LMPC/functions.jl")
include("barc_lib/LMPC/MPC_models.jl")
include("barc_lib/LMPC/coeffConstraintCost.jl")
include("barc_lib/LMPC/solveMpcProblem.jl")
include("barc_lib/simModel.jl")

# This function is called whenever a new state estimate is received.
# It saves this estimate in oldTraj and uses it in the MPC formulation (see in main)
function SE_callback(msg::pos_info,acc_f::Array{Float64},lapStatus::LapStatus,posInfo::PosInfo,mpcSol::MpcSol,oldTraj::OldTrajectory,trackCoeff::TrackCoeff,z_est::Array{Float64,1},x_est::Array{Float64,1},rhoEst::Array{Float64},epsiRef::Array{Float64})         # update current position and track data
    # update mpc initial condition
    v_abs                   =   sqrt(msg.v_x^2+msg.v_y^2)
    z_est[:]                =   [msg.s,msg.ey,msg.epsi,v_abs,rhoEst[1],epsiRef[1],acc_f[1]]             # use z_est as pointer
    x_est[:]                  = [msg.x,msg.y,msg.psi,msg.v]
    trackCoeff.coeffCurvature = msg.coeffCurvature
    # TODO: Check consequences of changing the order in z_est and in zCurr

    # check if lap needs to be switched
    if z_est[1] <= lapStatus.s_lapTrigger && lapStatus.switchLap
        oldTraj.idx_end[lapStatus.currentLap] = oldTraj.count[lapStatus.currentLap]
        oldTraj.oldCost[lapStatus.currentLap] = oldTraj.idx_end[lapStatus.currentLap] - oldTraj.idx_start[lapStatus.currentLap] #TODO: Fix costs
        lapStatus.currentLap += 1
        lapStatus.nextLap = true
        lapStatus.switchLap = false
    elseif z_est[1] > lapStatus.s_lapTrigger
        lapStatus.switchLap = true
    end

    # save current state in oldTraj
    oldTraj.oldTraj[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = z_est
    oldTraj.oldInput[oldTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = [msg.u_a,msg.u_df]
    oldTraj.oldTimes[oldTraj.count[lapStatus.currentLap],lapStatus.currentLap] = to_sec(msg.header.stamp)
    oldTraj.count[lapStatus.currentLap] += 1

end

function computeCost!(mpcTraj::MpcTrajectory,lapStatus::LapStatus,posInfo::PosInfo,mpcParams::MpcParams,modelParams::ModelParams) 

    lapNum = lapStatus.currentLap-1
    #get weights
    Q_term          = mpcParams.Q_term
    Q_term_cost     = mpcParams.Q_term_cost       
    R               = mpcParams.R                             
    QderivZ         = mpcParams.QderivZ          
    QderivU         = mpcParams.QderivU               
    Q_modelError    = mpcParams.Q_modelError
    rhoRef         = 1.0/modelParams.l_A
    # get trajectories (extract only relevant data)
    data_end = mpcTraj.idx_end[lapNum]
    stateHistory = mpcTraj.closedLoopSEY[1:data_end,:,lapNum]
    inputHistory = mpcTraj.inputHistory[1:data_end,:,lapNum]
    epsiRefHistory = mpcTraj.epsiRef[1:data_end,lapNum]
    # compute costs recursively, starting with the the last state with s < lapLength 
    # (note: for all other states that come after the finish line a cost of 0 is assumed (default value in array))
    for iii in data_end:-1:1
        s = stateHistory[iii,1]
        # -------------------------------------------------------------------------
        # determine stage cost components 
        if iii == 1
            derivCost = 0.0
            epsiRef = stateHistory[iii,3]

        else
            epsiRef = epsiRefHistory[iii-1]
            derivCost = 0.0
            for j = 1:5
                derivCost += QderivZ[j]*(stateHistory[iii,j] - stateHistory[iii-1,j])^2
            end
            for j = 1:3
                derivCost += QderivU[j]*(inputHistory[iii,j] - inputHistory[iii-1,j])^2
            end
        end

        controlCost = 0.0
        for j = 1:2
            controlCost += R[j]*(inputHistory[iii,j])^2
        end
        modelErrorCost = Q_modelError*(stateHistory[iii,3]-epsiRef)^2
        #modelErrorCost = Q_modelError*(stateHistory[iii,5]-rhoRef)^2
        # -------------------------------------------------------------------------

        if s >= posInfo.s_target
            currentStageCost = 0.0
        else
            currentStageCost = Q_term_cost + derivCost + controlCost + modelErrorCost  
        end 
        mpcTraj.cost[iii,2:end,lapNum] = [Q_term_cost;derivCost;controlCost;modelErrorCost;mpcTraj.cost[iii+1,1,lapNum]]

        # example: cost at t4 = stagecost at t4 + total cost of t5
        # total cost:
        mpcTraj.cost[iii,1,lapNum] =  mpcTraj.cost[iii+1,1,lapNum] + currentStageCost
    end

    return nothing

end
# this is only used to determine the chosen terminal states (for plotting)
function evaluateTerminalConstraint(mpcSol::MpcSol,mpcCoeff::MpcCoeff,mpcParams::MpcParams)
    N = mpcParams.N
    ParInt = mpcSol.ParInt
    ss = mpcCoeff.coeffConst #safe set approximation polynomial coefficients
    sF = mpcSol.z[N+1,1]

    # ey_F = (ParInt*(sF^5*ss[1,1,1]+sF^4*ss[2,1,1]+sF^3*ss[3,1,1]+sF^2*ss[4,1,1]+sF*ss[5,1,1]+ss[6,1,1]) + (1-ParInt)*(sF^5*ss[1,2,1]+sF^4*ss[2,2,1]+sF^3*ss[3,2,1]+sF^2*ss[4,2,1]+sF*ss[5,2,1]+ss[6,2,1]))

    # epsi_F = (ParInt*(sF^5*ss[1,1,2]+sF^4*ss[2,1,2]+sF^3*ss[3,1,2]+sF^2*ss[4,1,2]+sF*ss[5,1,2]+ss[6,1,2]) + (1-ParInt)*(sF^5*ss[1,2,2]+sF^4*ss[2,2,2]+sF^3*ss[3,2,2]+sF^2*ss[4,2,2]+sF*ss[5,2,2]+ss[6,2,2]))

    # v_F = (ParInt*(sF^5*ss[1,1,3]+sF^4*ss[2,1,3]+sF^3*ss[3,1,3]+sF^2*ss[4,1,3]+sF*ss[5,1,3]+ss[6,1,3]) + (1-ParInt)*(sF^5*ss[1,2,3]+sF^4*ss[2,2,3]+sF^3*ss[3,2,3]+sF^2*ss[4,2,3]+sF*ss[5,2,3]+ss[6,2,3]))

    # rho_F = (ParInt*(sF^5*ss[1,1,4]+sF^4*ss[2,1,4]+sF^3*ss[3,1,4]+sF^2*ss[4,1,4]+sF*ss[5,1,4]+ss[6,1,4]) + (1-ParInt)*(sF^5*ss[1,2,4]+sF^4*ss[2,2,4]+sF^3*ss[3,2,4]+sF^2*ss[4,2,4]+sF*ss[5,2,4]+ss[6,2,4]))
    xF = zeros(5)
    #xF = [sF;ey_F;epsi_F;v_F;rho_F]
    return xF
end
# This is the main function, it is called when the node is started.
function main()

    println("Starting LMPC node.")

    buffersize                  = 2000       # size of oldTraj buffers

    # Define and initialize variables
    # ---------------------------------------------------------------
    # General LMPC variables
    oldTraj                     = OldTrajectory() #m: data structure to save information about prev. trajectories
    posInfo                     = PosInfo() #m: information about current position s and s_target
    mpcCoeff                    = MpcCoeff() #m: coefficients for SS-traj and cost approximation
    lapStatus                   = LapStatus(1,1,false,false,0.3) #m: information about lap number and iteration as well as trigger and help variables
    mpcSol                      = MpcSol() #m: ax, df and other decision variables that come out of mpc
    trackCoeff                  = TrackCoeff()      # info about track (at current position, approximated)
    modelParams                 = ModelParams() #m: car dimensions, weight, actuator constraints, etc.
    mpcParams                   = MpcParams() #m: MPC horizon, weights, etc.
    mpcParams_pF                = MpcParams()       # for 1st lap (path following)
    mpcTraj                     = MpcTrajectory()   # data type to store the mpc trajectories and cost

    #m: Initialize all parameters with values defined in functions.jl
    InitializeParameters(mpcParams,mpcParams_pF,trackCoeff,modelParams,posInfo,oldTraj,mpcTraj,mpcCoeff,lapStatus,buffersize)
    # create two models: one for LMPC one for PF
    mdl    = MpcModel(mpcParams,mpcCoeff,modelParams,trackCoeff,posInfo)
    mdl_pF = MpcModel_pF(mpcParams_pF,modelParams,trackCoeff)

    max_N = max(mpcParams.N,mpcParams_pF.N)
    # ROS-specific variables
    z_est                       = zeros(7)          # this is a buffer that saves current state information (xDot, yDot, psiDot, ePsi, eY, s)
    x_est                       = zeros(4)          # this is a buffer that saves further state information (x, y, psi, v)
    coeffX                      = zeros(9)          # buffer for coeffX (only logging)
    coeffY                      = zeros(9)          # buffer for coeffY (only logging)
    cmd                         = ECU()             # command type
    coeffCurvature_update       = zeros(trackCoeff.nPolyCurvature+1)

    # Logging variables
    log_coeff_Cost              = NaN*ones(mpcCoeff.order+1,2,10000)
    log_coeff_Const             = NaN*ones(mpcCoeff.order+1,2,4,400,30)
    log_sol_z                   = NaN*ones(max_N+1,7,10000)
    log_sol_u                   = NaN*ones(max_N,3,10000)
    log_curv                    = zeros(10000,trackCoeff.nPolyCurvature+1)
    log_state_x                 = zeros(10000,4)
    log_coeffX                  = zeros(10000,9)
    log_coeffY                  = zeros(10000,9)
    log_t                       = zeros(10000,1)
    log_state                   = zeros(10000,7)
    log_cost                    = zeros(10000,8)
    log_cmd                     = zeros(10000,3)
    log_step_diff               = zeros(10000,5)
    log_t_solv                  = zeros(10000)
    log_numIter                 = zeros(Int64,30)
    log_IterSwitch              = zeros(Int64,30)
    log_sol_status              = Array(Symbol,10000)
    log_mpc_sol_z               = zeros(400,max_N+1,7,30)
    acc_f                       = [0.0]
    rhoRef                      = 1.0/modelParams.l_A
    rhoEst                      = [rhoRef]                # for pathfollowing, will be later overwriten in mpc laps
    epsiRef                     = [0.0]
    # Initialize ROS node and topics
    init_node("mpc_traj")
    loop_rate = Rate(1/modelParams.dt)
    pub = Publisher("ecu", ECU, queue_size=1)::RobotOS.Publisher{barc.msg.ECU}
    # The subscriber passes arguments (coeffCurvature and z_est) which are updated by the callback function:
    s1 = Subscriber("pos_info", pos_info, SE_callback, (acc_f,lapStatus,posInfo,mpcSol,oldTraj,trackCoeff,z_est,x_est,rhoEst,epsiRef),queue_size=50)::RobotOS.Subscriber{barc.msg.pos_info}
    # Note: Choose queue size long enough so that no pos_info packets get lost! They are important for system ID!

    run_id = get_param("run_id")
    println("Finished initialization.")
    
    # buffer in current lap
    zCurr                       = zeros(10000,7)    # contains state information in w lap (max. 10'000 steps)
    uCurr                       = zeros(10000,3)    # contains input information
    step_diff                   = zeros(5)

    # Specific initializations:
    lapStatus.currentLap    = 1
    lapStatus.currentIt     = 1
    #posInfo.s_target        = 19.11                    
    k                       = 0                       # overall counter for logging
    
    mpcSol.z = zeros(11,4)
    mpcSol.u = zeros(10,2)
    mpcSol.a_x = 0
    mpcSol.d_f = 0
    
    # Precompile coeffConstraintCost:
    xfRange = zeros(2,2)
    selected_Laps = zeros(Int64,2)
    mpcTraj.closedLoopSEY[1:400,1,1] = linspace(0,posInfo.s_target,400)
    mpcTraj.closedLoopSEY[1:400,1,2] = linspace(0,posInfo.s_target,400)
    mpcTraj.count[1:2] = 400
    posInfo.s = posInfo.s_target/2
    lapStatus.currentLap = 3
    (xfRange,selected_Laps) = coeffConstraintCost(mpcTraj,mpcCoeff,posInfo,mpcParams,lapStatus)
    lapStatus.currentLap = 1
    mpcTraj.closedLoopSEY[1:400,1,1] = zeros(400,1)
    mpcTraj.closedLoopSEY[1:400,1,2] = zeros(400,1)
    mpcTraj.count[1:2] = 1

    posInfo.s = 0
    println("Precompiling of coeffConstraintCost() completed!")

    #TODO: Precompile computeCost! ?

    uPrev = zeros(10,2)     # saves the last 10 inputs (1 being the most recent one)

    n_pf = 2               # number of first path-following laps (needs to be at least 2)

    acc0 = 0.0
    opt_count = 0
    counter = 0
    # Start node
    while ! is_shutdown()
        if z_est[1] > 0         # check if data has been received (s > 0)
            # ============================= PUBLISH COMMANDS =============================
            # This is done at the beginning of the lap because this makes sure that the command is published 0.1s after the state has been received
            # This guarantees a constant publishing frequency of 10 Hz
            # (The state can be predicted by 0.1s)
            cmd.header.stamp = get_rostime()
            # cmd.motor = convert(Float32,mpcSol.a_x)
            # cmd.servo = convert(Float32,mpcSol.d_f)
            publish(pub, cmd)
            # ============================= Initialize iteration parameters =============================
            i                           = lapStatus.currentIt           # current iteration number, just to make notation shorter
            zCurr[i,:]                  = copy(z_est)                   # update state information
            posInfo.s                   = zCurr[i,1]                    # update position info
            #trackCoeff.coeffCurvature   = copy(coeffCurvature_update)

            # ============================= Pre-Logging (before solving) ================================
            log_t[k+1]                  = to_sec(get_rostime())         # time is measured *before* solving (more consistent that way)
            # if size(mpcSol.z,2) == 4                                    # find 1-step-error
            #     step_diff = ([mpcSol.z[2,4], 0, 0, mpcSol.z[2,3], mpcSol.z[2,2]]-[norm(zCurr[i,1:2]), 0, 0, zCurr[i,4], zCurr[i,5]])
            # else
            #     step_diff = (mpcSol.z[2,1:5][:]-zCurr[i,1:5][:])
            # end
            # log_step_diff[k+1,:]          = step_diff

            # ======================================= Lap trigger =======================================
            if lapStatus.nextLap                # if we are switching to the next lap...
                println("Finishing one lap at iteration ",i)
                # Important: lapStatus.currentIt is now the number of points up to s > s_target -> -1 in saveOldTraj
                zCurr[1,:] = zCurr[i,:]         # copy current state
                mpcTraj.idx_end[lapStatus.currentLap-1] = mpcTraj.count[lapStatus.currentLap-1] - 1  #save the number of mpc steps per lap
                computeCost!(mpcTraj,lapStatus,posInfo,mpcParams,modelParams)           # compute costs for states in the last lap recursively
                log_IterSwitch[lapStatus.currentLap-1] = k-1

                i                     = 1
                lapStatus.currentIt   = 1       # reset current iteration
                lapStatus.nextLap = false

                # Set warm start for new solution (because s shifted by s_target)
                if lapStatus.currentLap <= n_pf
                    setvalue(mdl_pF.z_Ol[:,1],mpcSol.z[:,1]-posInfo.s_target)
                elseif lapStatus.currentLap == n_pf+1
                    setvalue(mdl.z_Ol[:,1],mpcSol.z[1:mpcParams.N+1,1]-posInfo.s_target)
                    setvalue(mdl.z_Ol[:,2],mpcSol.z[1:mpcParams.N+1,2])
                    setvalue(mdl.z_Ol[:,3],mpcSol.z[1:mpcParams.N+1,3])
                    setvalue(mdl.z_Ol[:,4],mpcSol.z[1:mpcParams.N+1,4])
                    setvalue(mdl.z_Ol[:,5],mpcSol.z[1:mpcParams.N+1,5])
                    #m: setvalue(mdl.u_Ol,mpcSol.u[1:mpcParams.N,:])
                elseif lapStatus.currentLap > n_pf+1
                    setvalue(mdl.z_Ol[:,1],mpcSol.z[:,1]-posInfo.s_target)
                end
            end

            #  ======================================= Calculate input =======================================
            #println("*** NEW ITERATION # ",i," ***")
            println("Current Lap: ", lapStatus.currentLap, ", It: ", lapStatus.currentIt)
            #println("State Nr. ", i, "    = ", z_est)
            #println("s               = $(posInfo.s)")
            #println("s_total         = $(posInfo.s%posInfo.s_target)")

            # Find coefficients for cost and constraints
            if lapStatus.currentLap > n_pf
                tic()
                (xfRange,selected_Laps) = coeffConstraintCost(mpcTraj,mpcCoeff,posInfo,mpcParams,lapStatus)
                mpcTraj.xfRange[mpcTraj.count[lapStatus.currentLap],:,:,lapStatus.currentLap] = xfRange
                mpcTraj.selected_Laps[mpcTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = selected_Laps
                tt = toq()
            end

            #println("Starting solving.")
            # Solve the MPC problem
            tic()
            if lapStatus.currentLap <= n_pf
                z_pf = [zCurr[i,1:4]';acc0]        # use kinematic model and its states
                solveMpcProblem_pathFollow(mdl_pF,mpcSol,mpcParams_pF,trackCoeff,posInfo,modelParams,z_pf,uPrev)
                acc_f[1] = mpcSol.z[1,5]
                acc0 = mpcSol.z[2,5]
                epsiRef[1] = mpcSol.z[2,3] #same as epsi
                z0 = [mpcSol.z[1,1:4]';rhoRef;mpcSol.z[1,3];acc_f[1]]
                u0 = [mpcSol.u[1,1:2]';0.0]
            else                        # otherwise: use adaptive kinematic model
                zCurr[i,7] = acc0
                zLMPC = [zCurr[i,1:6]';acc0]    
                solveMpcProblem(mdl,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zLMPC,uPrev)
                acc0 = mpcSol.z[2,7]
                acc_f[1] = mpcSol.z[1,7]
                rhoEst[1] = mpcSol.z[2,5]
                epsiRef[1] = mpcSol.z[2,6]
                z0 = mpcSol.z[1,:]
                u0 = mpcSol.u[1,:]
                mpcTraj.xfStates[mpcTraj.count[lapStatus.currentLap],:,lapStatus.currentLap] = evaluateTerminalConstraint(mpcSol,mpcCoeff,mpcParams)
                mpcTraj.epsiRef[mpcTraj.count[lapStatus.currentLap],lapStatus.currentLap] = mpcSol.z[2,6]
            end

            #  ======================================= Save States for Safe Set =======================================
            # save states that were used in the mpc for safe set computation
            # if necessary: append to end of previous lap
            mpcTraj.closedLoopSEY[mpcTraj.count[lapStatus.currentLap],1:7,lapStatus.currentLap] = z0
            mpcTraj.inputHistory[mpcTraj.count[lapStatus.currentLap],1:3,lapStatus.currentLap] = u0
            mpcTraj.count[lapStatus.currentLap] += 1
            # also add some states to the end of previous lap (for safe set approximation)
            if lapStatus.currentLap > 1 && z0[1] < 8.0
                mpcTraj.closedLoopSEY[mpcTraj.count[lapStatus.currentLap-1],1:7,lapStatus.currentLap-1] = z0
                mpcTraj.closedLoopSEY[mpcTraj.count[lapStatus.currentLap-1],1,lapStatus.currentLap-1] += posInfo.s_target
                mpcTraj.count[lapStatus.currentLap-1] += 1
            end




            log_t_solv[k+1] = toq()


            cmd.motor = convert(Float32,mpcSol.a_x)
            cmd.servo = convert(Float32,mpcSol.d_f)    

            # Write current input information
            uCurr[i,:] = [mpcSol.a_x mpcSol.d_f mpcSol.phi]
            zCurr[i,1] = posInfo.s%posInfo.s_target   # save absolute position in s (for oldTrajectory)

            uPrev = circshift(uPrev,1) #m: shift data by one up
            uPrev[1,:] = uCurr[i,1:2] #m: change the first entry 
            #println("Finished solving, status: $(mpcSol.solverStatus), u = $(uCurr[i,:]), t = $(log_t_solv[k+1]) s")

            # Logging
            # ---------------------------
            k = k + 1       # counter
            log_sol_status[k]       = mpcSol.solverStatus
            log_state[k,:]          = zCurr[i,:]
            log_cmd[k+1,:]          = uCurr[i,:]                    # the command is going to be pubished in the next iteration
            log_coeff_Cost[:,:,k]   = mpcCoeff.coeffCost
            log_coeff_Const[:,:,:,lapStatus.currentIt,lapStatus.currentLap] = mpcCoeff.coeffConst
            log_cost[k,:]           = mpcSol.cost
            log_curv[k,:]           = trackCoeff.coeffCurvature
            log_state_x[k,:]        = x_est
            if size(mpcSol.z,2) == 5 #Path Following Case
                log_sol_z[1:mpcParams_pF.N+1,1:4,k]   = mpcSol.z[:,1:4]        # log s, ey, epsi, v
                log_sol_z[1:mpcParams_pF.N+1,5,k]     = rhoEst[1]*ones(mpcParams_pF.N+1,1)        # rhoEst, epsiRef
                log_sol_z[1:mpcParams_pF.N+1,6,k]     = mpcSol.z[:,3]        # in path following case there is no difference between epsi and epsiRef
                log_sol_z[1:mpcParams_pF.N+1,7,k]     = mpcSol.z[:,5]        # in path following filter state is number 5

                log_sol_u[1:mpcParams_pF.N,1:2,k]         = mpcSol.u        # log (a, dF)
            else #LMPC case
                log_mpc_sol_z[lapStatus.currentIt,1:mpcParams.N+1,1:7,lapStatus.currentLap] = mpcSol.z
                log_sol_z[1:mpcParams.N+1,1:7,k]        = mpcSol.z
                log_sol_u[1:mpcParams.N,:,k]            = mpcSol.u
            end

            # Count one up:
            lapStatus.currentIt += 1
            log_numIter[lapStatus.currentLap] = counter
            counter +=1
        else
            println("No estimation data received!")
        end
        rossleep(loop_rate)
    end
    # Save simulation data to file

    log_path = "$(homedir())/simulations/output-LMPC-$(run_id[1:4]).jld"
    if isfile(log_path)
        log_path = "$(homedir())/simulations/output-LMPC-$(run_id[1:4])-2.jld"
        warn("Warning: File already exists.")
    end
    save(log_path,"oldTraj",oldTraj,"state",log_state[1:k,:],"t",log_t[1:k],"sol_z", log_sol_z,"sol_u",log_sol_u[:,:,1:k],
        "cost",log_cost[1:k,:],"curv",log_curv[1:k,:],"coeffCost",log_coeff_Cost,"coeffConst",log_coeff_Const,
        "x_est",log_state_x[1:k,:],"cmd",log_cmd[1:k,:],"t_solv",log_t_solv[1:k],"sol_status",log_sol_status[1:k],"numIter",log_numIter,"mpcTraj",mpcTraj,"mpcParams",mpcParams,"mpcParams_pF",mpcParams_pF,"log_mpc_sol_z",log_mpc_sol_z,"IterSwitch",log_IterSwitch)
    println("Exiting LMPC node. Saved data to $log_path.")

end

if ! isinteractive()
    main()   
end



# Sequence within one iteration:
# 1. Publish commands from last iteration (because the car is in real *now* where we thought it was before (predicted state))
# 2. Receive new state information
# 3. Check if we've crossed the finish line and if we have, switch lap number and save old trajectories
# 4. (If in 3rd lap): Calculate coefficients
# 5. Calculate MPC commands (depending on lap) which are going to be published in the next iteration
# 6. (Do some logging)


# Definitions of variables:
# zCurr contains all state information from the beginning of the lap (first s >= 0) to the current state i
# uCurr -> same as zCurr for inputs
# generally: zCurr[i+1] = f(zCurr[i],uCurr[i])

# zCurr[1] = v_x
# zCurr[2] = v_y
# zCurr[3] = psiDot
# zCurr[4] = ePsi
# zCurr[5] = eY
# zCurr[6] = s
