module debugModule


    using JLD
    using PyPlot
    using PyCall
      @pyimport matplotlib.animation as animation
      @pyimport matplotlib as mpl
      @pyimport matplotlib.patches as patches
      using JLD, ProfileView
    export Measurements, plotStates, animateState2
    # pos_info[1]  = s
    # pos_info[2]  = eY
    # pos_info[3]  = ePsi
    # pos_info[4]  = v
    # pos_info[5]  = s_start
    # pos_info[6]  = x
    # pos_info[7]  = y
    # pos_info[8]  = v_x
    # pos_info[9]  = v_y
    # pos_info[10] = psi
    # pos_info[11] = psiDot
    # pos_info[12] = x_raw
    # pos_info[13] = y_raw
    # pos_info[14] = psi_raw
    # pos_info[15] = v_raw
    # pos_info[16] = psi_drift
    # pos_info[17] = a_x
    # pos_info[18] = a_y

    include("../workspace/src/barc/src/barc_lib/classes.jl")

    type Measurements{T}
        i::Int64                # measurement counter
        t::Array{Float64}       # time data (when it was received by this recorder)
        t_msg::Array{Float64}   # time that the message was sent
        z::Array{T}             # measurement values
    end


    function plotStates(code::AbstractString,stateNum::Int64,selectedLaps::Array{Int64}=[0],plotOverK::Bool=true)
        # stateNum:

   
        PyPlot.close(222)

        log_path_sim = "$(homedir())/simulations/output-SIM-$(code).jld"
        log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
        log_path_LMPC   = "$(homedir())/simulations/output-LMPC-$(code).jld"

        d_rec       = load(log_path_record)
        d_lmpc      = load(log_path_LMPC)


        # load the data

        # from recorder node
        pos_info    = d_rec["pos_info"]

        # from lmpc node
        oldTraj     = d_lmpc["oldTraj"]
        t           = d_lmpc["t"]
        state       = d_lmpc["state"]
        sol_z       = d_lmpc["sol_z"]
        sol_u       = d_lmpc["sol_u"]
        cost        = d_lmpc["cost"]
        numIter     = d_lmpc["numIter"]
        # get t0
        #t0 = pos_info.t[1]

        # design settings
        fsize = 30


        # determine how many laps were completed in the file
        numLaps = length(find(x -> x>1, oldTraj.oldCost))
        if numLaps < 1; error("At least one lap has to be completed succesfully!") end

        # handle case where user wants to plot all laps, if so add all lap numbers to array
        if selectedLaps == [0] 
            selectedLaps = 1:1:numLaps
        end 


        # define new figure
        fig = PyPlot.figure(222,facecolor="white")
        PyPlot.hold(true)#
        #(nothing, axarr) = subplots(2, sharex=true)
        #ax1 = axarr[1]
        #ax2 = axarr[2]
        ax1 = fig[:add_subplot](2, 1, 1)
        ax2 = fig[:add_subplot](2, 1, 2)
        fig[:subplots_adjust](hspace=.5)
        # iterate through laps, select data and plot
        for iii in selectedLaps

            # determine data range for measured data
            start_meas = oldTraj.idx_start[iii]
            if start_meas == 0; start_meas = 1 end
            end_meas = oldTraj.idx_end[iii]

            # determine data range for mpc data
            if iii == 1
                start_mpc = 1
            else
                start_mpc = numIter[iii-1] + 2 
            end
            end_mpc = numIter[iii] + 1
            # plot over k or over s
            if plotOverK   
                stateData_meas = oldTraj.oldTraj[start_meas:end_meas,stateNum,iii][:]
                xdata_meas = 1:1:length(stateData_meas)

                stateData_mpc = sol_z[1,stateNum,start_mpc:end_mpc][:]
                xdata_mpc = 1:1:length(stateData_mpc)

            else
                stateData_meas = oldTraj.oldTraj[start_meas:end_meas,stateNum,iii][:]
                xdata_meas = oldTraj.oldTraj[start_meas:end_meas,1,iii][:]

                stateData_mpc = sol_z[1,stateNum,start_mpc:end_mpc][:]
                xdata_mpc = sol_z[1,1,start_mpc:end_mpc][:]
            end
            ax1[:plot](xdata_meas,stateData_meas, label="Lap: $(iii)", linewidth=1.0)
            ax2[:plot](xdata_mpc,stateData_mpc, label="Lap: $(iii)", linewidth=1.0)
        end
        if plotOverK
            ax1[:set_xlabel](L"k",fontsize=fsize)
            ax2[:set_xlabel](L"k",fontsize=fsize)
        else
            ax1[:set_xlabel](L"s",fontsize=fsize)
            ax2[:set_xlabel](L"s",fontsize=fsize)        
        end

        ylabels = ["s","ey","epsi","v","rhoEst","epsiRef","filter"]
        ax1[:set_ylabel](ylabels[stateNum]*" (meas)",fontsize=fsize)
        ax2[:set_ylabel](ylabels[stateNum]*" (mpc)",fontsize=fsize)
        ax1[:grid]()
        ax2[:grid]()
        ax1[:legend]()
        ax2[:legend]()
        PyPlot.hold(false)
    end


    function eval_sim(code::AbstractString)
        log_path_sim = "$(homedir())/simulations/output-SIM-$(code).jld"
        log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
        d_sim = load(log_path_sim)
        d_rec = load(log_path_record)

        imu_meas    = d_sim["imu_meas"]
        gps_meas    = d_sim["gps_meas"]
        z           = d_sim["z"]
        cmd_log     = d_sim["cmd_log"]
        slip_a      = d_sim["slip_a"]
        pos_info    = d_rec["pos_info"]
        vel_est     = d_rec["vel_est"]

        t0 = pos_info.t[1]
        track = create_track(0.4)

        figure()
        ax1=subplot(311)
        plot(z.t-t0,z.z,"-*")
        title("Real states")
        grid()
        legend(["x","y","v_x","v_y","psi","psi_dot","a","d_f"])
        subplot(312,sharex=ax1)
        plot(cmd_log.t-t0,cmd_log.z,"-*")
        title("Inputs")
        grid()
        legend(["u","d_f"])
        subplot(313,sharex=ax1)
        plot(slip_a.t-t0,slip_a.z,"-*")
        title("Slip angles")
        grid()
        legend(["a_f","a_r"])

        figure()
        plot(z.z[:,1],z.z[:,2],"-",gps_meas.z[:,1],gps_meas.z[:,2],".",pos_info.z[:,6],pos_info.z[:,7],"-")
        plot(track[:,1],track[:,2],"b.",track[:,3],track[:,4],"r-",track[:,5],track[:,6],"r-")
        grid(1)
        title("x-y-view")
        axis("equal")
        legend(["Real state","GPS meas","Estimate"])
        
        figure()
        title("Comparison of psi")
        plot(imu_meas.t-t0,imu_meas.z,"-x",z.t-t0,z.z[:,5:6],pos_info.t-t0,pos_info.z[:,10:11],"-*")
        legend(["imu_psi","imu_psi_dot","real_psi","real_psi_dot","est_psi","est_psi_dot"])
        grid()

        figure()
        title("Comparison of v")
        plot(z.t-t0,z.z[:,3:4],z.t-t0,sqrt(z.z[:,3].^2+z.z[:,4].^2),pos_info.t-t0,pos_info.z[:,8:9],"-*",pos_info.t-t0,sqrt(pos_info.z[:,8].^2+pos_info.z[:,9].^2),"-*",vel_est.t-t0,vel_est.z)
        legend(["real_xDot","real_yDot","real_v","est_xDot","est_yDot","est_v","v_x_meas"])
        grid()

        figure()
        title("Comparison of x,y")
        plot(z.t-t0,z.z[:,1:2],pos_info.t-t0,pos_info.z[:,6:7],"-*",gps_meas.t-t0,gps_meas.z)
        legend(["real_x","real_y","est_x","est_y","meas_x","meas_x"])
        grid()
    end

     function animateState2(code::AbstractString,iii::Int64,stateNum::Int64=2,plotOverK::Bool=true,offset::Int64=0,speed::Int64=500 )
    if iii < 3; error("Please specify Iteration from LMPC with iterNum > 2") end
    if offset < 0; error("Offset must be a positive integer") end
    if !(2 <= stateNum <= 5); error("Specify a state number between 2 and 6") end
    
    log_path_sim = "$(homedir())/simulations/output-SIM-$(code).jld"
    log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
    log_path_LMPC   = "$(homedir())/simulations/output-LMPC-$(code).jld"

    d_rec       = load(log_path_record)
    d_lmpc      = load(log_path_LMPC)


    # load the data

    # from recorder node
    pos_info    = d_rec["pos_info"]

    # from lmpc node
    oldTraj     = d_lmpc["oldTraj"]
    mpcTraj     = d_lmpc["mpcTraj"]
    t           = d_lmpc["t"]
    state       = d_lmpc["state"]
    sol_z       = d_lmpc["sol_z"]
    sol_u       = d_lmpc["sol_u"]
    cost        = d_lmpc["cost"]
    numIter     = d_lmpc["numIter"]



    # TODO: Remove hard code
    N_PF = 15
    N = 12 

    # Define array of state names
    stateNames = ["s","ey","epsi","v","rho","epsiRef","filter"]
    pygui(true)
    fig = figure(facecolor="white",figsize=(15,10))
    PyPlot.grid(true)
    PyPlot.ylabel(stateNames[stateNum],fontsize=20)
    PyPlot.title("Trajectory and open loop for $(stateNames[stateNum])",fontsize=20)

    steps_prev = length(mpcTraj.closedLoopSEY[:,1,iii-1])
    steps_cl = length(mpcTraj.closedLoopSEY[:,1,iii])
    # Define lines
    state_prev = plot([],[],"g-",linewidth=2.0,label="PF C.l.")[1]
    state_closedloop = plot([],[],"k-",linewidth=2.0,label="Sys C.L.")[1]
    #state_mismatch = plot([],[],"c--",linewidth=2.0,label="MPC C.L.")[1]
    state_openloop = plot([],[],"b-+",linewidth=2.0,markersize=10.0, label="O.L. Pred")[1]
    state_approx = plot([],[],"r-",linewidth=2.0, label="SS approx")[1]
    state_approxPoints = plot([],[],"r+",markersize=10.0)[1]
    PyPlot.legend(fontsize=20,bbox_to_anchor=(1.1, 1.1))

    kMin = 1
    kMax = steps_prev
    sMin = 0
    sMax = maximum(mpcTraj.closedLoopSEY[:,1,iii])
    eyMin = - 0.6
    eyMax = + 0.6
    epsiMin = minimum(mpcTraj.closedLoopSEY[:,3,iii]) + 0.3
    epsiMax = maximum(mpcTraj.closedLoopSEY[:,3,iii]) + 0.3
    vMin = - 0.4
    vMax = 3.4
    rhoMin = -0.1
    rhoMax = maximum(mpcTraj.closedLoopSEY[:,5,iii]) + 5
    

    ylimBtm = [sMin,eyMin,epsiMin,vMin,rhoMin,epsiMin,-10]
    ylimTop = [sMax,eyMax,epsiMax,vMax,rhoMax,epsiMax,10]
    axes =gca()

    myMin = ylimBtm[stateNum]
    myMax = ylimTop[stateNum]
    axes[:set_ylim]([myMin,myMax])

    if plotOverK
      axes[:set_xlim]([kMin,kMax])
      PyPlot.xlabel("k",fontsize=20)
    else
      axes[:set_xlim]([sMin,sMax])
      PyPlot.xlabel("s",fontsize=20)
    end


    function init_plot()
        # plot background: previous trajectory (SS) and closed loop trajectory
        if plotOverK
          state_prev[:set_data]([1:1:steps_prev],[mpcTraj.closedLoopSEY[:,stateNum,iii-1]])
          state_closedloop[:set_data]([1:1:steps_cl],[mpcTraj.closedLoopSEY[:,stateNum,iii]])
          #state_mismatch[:set_data]([1:1:steps_cl],[Iterations[iii].openLoopSEY[stateNum,2,:]]) #FIXME Also check here if ytrue time index used or it has to be t-1
        else
          state_prev[:set_data]([Iterations[iii-1].closedLoopSEY[1,:]],[mpcTraj.closedLoopSEY[:,stateNum,iii-1]])
          state_closedloop[:set_data]([mpcTraj.closedLoopSEY[:,1,iii]],[mpcTraj.closedLoopSEY[:,stateNum,iii]])
          #state_mismatch[:set_data]([Iterations[iii].openLoopSEY[1,2,:]],[Iterations[iii].openLoopSEY[stateNum,2,:]])
        end

      return (state_prev,state_closedloop,nothing)
    end

    function animate(frameNum)
      frameNum += 1
      k = frameNum + offset

      if plotOverK
        # Update open loop trajectory data
        state_openloop[:set_data]([k:1:k+N],[Iterations[iii].openLoopSEY[stateNum,:,k]])

        # Update s-ey data points of safe set trajectory
        xfRangeMin = Int64(Iterations[iii].xfRange[1,k])
        xfRangeMax = Int64(Iterations[iii].xfRange[2,k])
        sRange = Iterations[iii-1].openLoopSEY[1,xfRangeMin:1:xfRangeMax]
        stateRange = Iterations[iii-1].openLoopSEY[stateNum,xfRangeMin:1:xfRangeMax]
        state_approxPoints[:set_data]([xfRangeMin:1:xfRangeMax],[stateRange])

        # evaluate polynomial approximation of safe set points and update data
        ss = Iterations[iii].xfApproxCoeff[:,stateNum-1,k]
        statePoly = sRange.^5 *ss[1] + sRange.^4 *ss[2] + sRange.^3 *ss[3] + sRange.^2 *ss[4] + sRange*ss[5] + ss[6]
        state_approx[:set_data]([xfRangeMin:1:xfRangeMax],[statePoly])
      else
        # Update open loop trajectory data
        state_openloop[:set_data]([Iterations[iii].openLoopSEY[1,:,k]],[Iterations[iii].openLoopSEY[stateNum,:,k]])

        # Update s-ey data points of safe set trajectory
        xfRangeMin = Int64(Iterations[iii].xfRange[1,k])
        xfRangeMax = Int64(Iterations[iii].xfRange[2,k])
        sRange = Iterations[iii-1].openLoopSEY[1,xfRangeMin:1:xfRangeMax]
        stateRange = Iterations[iii-1].openLoopSEY[stateNum,xfRangeMin:1:xfRangeMax]
        state_approxPoints[:set_data]([sRange],[stateRange])

        # evaluate polynomial approximation of safe set points and update data
        ss = Iterations[iii].xfApproxCoeff[:,stateNum-1,k]
        statePoly = sRange.^5 *ss[1] + sRange.^4 *ss[2] + sRange.^3 *ss[3] + sRange.^2 *ss[4] + sRange*ss[5] + ss[6]
        state_approx[:set_data]([sRange],[statePoly])
      end
      return (state_openloop,state_approx,state_approxPoints,nothing)
    end
    NumFrames = length(Iterations[iii].openLoopSEY[1,:]) - offset
    anim = animation.FuncAnimation(fig, animate, init_func=init_plot, frames=NumFrames, interval=speed,repeat=false)
  end
end #MODULE
    # add_curve(theta,30,0)
    # add_curve(theta,20,-pi/6)