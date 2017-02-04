include("ColorModule.jl")

module debugModule

using JLD
using PyPlot
using PyCall
using ColorModule
@pyimport matplotlib.animation as animation
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
export Measurements, plotStates, animateState2, create_track
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

    # Load colors from user defined object Color Manager
    CM = ColorManager()

    # define new figure
    fig = PyPlot.figure(facecolor="white")
    PyPlot.hold(true)#
    #(nothing, axarr) = subplots(2, sharex=true)
    #ax1 = axarr[1]
    #ax2 = axarr[2]
    # ax1 = fig[:add_subplot](2, 1, 1)
    # ax2 = fig[:add_subplot](2, 1, 2)
    # fig[:subplots_adjust](hspace=.5)
    # iterate through laps, select data and plot
    for iii in selectedLaps

      # # determine data range for measured data
      # start_meas = oldTraj.idx_start[iii]
      # if start_meas == 0; start_meas = 1 end
      # end_meas = oldTraj.idx_end[iii]

      # determine data range for mpc data
      if iii == 1
        start_mpc = 1
      else
        start_mpc = numIter[iii-1] + 2 
      end
    end_mpc = numIter[iii] + 1
    # plot over k or over s
    if plotOverK   
      # stateData_meas = oldTraj.oldTraj[start_meas:end_meas,stateNum,iii][:]
      # xdata_meas = 1:1:length(stateData_meas)

      stateData_mpc = sol_z[1,stateNum,start_mpc:end_mpc][:]
      xdata_mpc = 1:1:length(stateData_mpc)

    else
      # stateData_meas = oldTraj.oldTraj[start_meas:end_meas,stateNum,iii][:]
      # xdata_meas = oldTraj.oldTraj[start_meas:end_meas,1,iii][:]

      stateData_mpc = sol_z[1,stateNum,start_mpc:end_mpc][:]
      xdata_mpc = sol_z[1,1,start_mpc:end_mpc][:]
    end
    #plot(xdata_meas,stateData_meas, label="Lap: $(iii)", linewidth=1.0)
    plot(xdata_mpc,stateData_mpc, label="Lap: $(iii)", linewidth=2.0,color=getColor(CM))
  end
  ax = gca()
  if plotOverK
    ax[:set_xlabel](L"k",fontsize=fsize)
  else
    ax[:set_xlabel](L"s",fontsize=fsize)
  end

  ylabels = ["s","ey","epsi","v","rhoEst","epsiRef","filter"]
  ax[:set_ylabel](ylabels[stateNum]*" (mpc)",fontsize=fsize)
  ax[:grid]()
  # ax2[:grid]()
  ax[:legend]()
  # ax2[:legend]()
  PyPlot.hold(false)
end

# plots the measured  x-y position of the car over the racetrack (lap wise)
function plotPositions(code::AbstractString)

  # create racetrack plot
  (xt, yt, xt_l, yt_l, xt_r, yt_r) = create_track(0.5)
  fig = figure(999,facecolor="white")
  plot(xt,yt,"k--")
  hold(true)
  plot(xt_l,yt_l,"k",xt_r,yt_r,"k",linewidth=2.0)
  grid(true)

  # load data from JLD files
  log_path_LMPC   = "$(homedir())/simulations/output-LMPC-$(code).jld"
  log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
  d_lmpc      = load(log_path_LMPC)
  d_rec = load(log_path_record)

  # determine how many laps were completed in the file
  state_x    = d_lmpc["x_est"]
  numIter    = d_lmpc["numIter"]
  IterSwitch = d_lmpc["IterSwitch"]

  CM = ColorManager()
  # determine how many laps were completed in the file
  numLaps = length(find(x -> x>1, numIter))
  if numLaps < 1; error("At least one lap has to be completed succesfully!") end

    for iii = 1:numLaps
      if iii == 1
        datarange = 1:IterSwitch[iii]
      else
        datarange = IterSwitch[iii-1]+1:IterSwitch[iii]
      end
      xpos = state_x[datarange,1]
      ypos = state_x[datarange,2]
      plot(xpos,ypos,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
    end
    xlabel(L"x")
    ylabel(L"y")
    legend()
    ax = gca()
    ax[:set_aspect]("equal")

  end

  function plotCost(code::AbstractString)

    # load data from JLD files
    log_path_LMPC     = "$(homedir())/simulations/output-LMPC-$(code).jld"
    d_lmpc          = load(log_path_LMPC)

    # determine how many laps were completed in the file
    mpcTraj    = d_lmpc["mpcTraj"]
    numIter    = d_lmpc["numIter"]
    costMPC    = d_lmpc["cost"]


    # determine how many laps were completed in the file
    numLaps = length(find(x -> x>1, numIter))
    if numLaps < 1; error("At least one lap has to be completed succesfully!") end

      CM = ColorManager()
      #toggle_mum(CM)
      # plot data  
      fig = figure(888,facecolor="white")
      title("Cost for each lap",fontsize=30)
      for iii = 1:numLaps
        sdata = mpcTraj.closedLoopSEY[1:mpcTraj.idx_end[iii],1,iii]
        cost = mpcTraj.cost[1:mpcTraj.idx_end[iii],1,iii]
        plot(sdata,cost,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
        hold(true)
      end
      grid(true)
      xlabel(L"s",fontsize=30)
      ylabel(L"cost",fontsize=30)
      legend()
      resetColorCycle(CM)
      # second plot over k
      fig2 = figure(444,facecolor="white")
      # plot the cost that were calculated after one lap was completed
      subplot(2,1,1)
      legendEntries = ["Total","Qterm_cost","derivCost","controlCost","modelErrorCost"]
      for jjj = 2:5
        cost = zeros(10000)
        startPoint = 1
        for iii = 1:numLaps
          dataLength = mpcTraj.idx_end[iii]
          #startPoint:1:dataLength
          cost[startPoint:1:startPoint+dataLength-1] = mpcTraj.cost[1:mpcTraj.idx_end[iii],jjj,iii] 
          startPoint = startPoint + dataLength
        end
        plot(cost[1:startPoint],label=legendEntries[jjj],color=getColor(CM))
        hold(true)
      end

      ax = gca()
      ax[:set_title]("Components of post determined cost")
      ax[:grid]()
      ax[:legend]()
      xlabel("k")
      ylabel("Cost component")


      resetColorCycle(CM)

      # show actual cost of mpc
      subplot(2,1,2)
      legendEntries = ["Total","stageCost","derivCost","controlCost","modelErrorCost","costZTerm","constZTerm","laneCost"]
      for jjj = 1:8
        plot(costMPC[:,jjj],label=legendEntries[jjj],color=getColor(CM))
        hold(true)
      end
      ax = gca()
      ax[:set_title]("Cost components MPC")
      ax[:grid]()
      ax[:legend]()
      xlabel("k")
      ylabel("Cost component")
    end

    function plotPredictionError(code::AbstractString)
    # load data from JLD files
      log_path_LMPC     = "$(homedir())/simulations/output-LMPC-$(code).jld"
      d_lmpc          = load(log_path_LMPC)

      # determine how many laps were completed in the file
      mpcTraj    = d_lmpc["mpcTraj"]
      numIter    = d_lmpc["numIter"]


      numLaps = length(find(x -> x>1, numIter))

      CM = ColorManager()
      # plot data  
      fig = figure(33,facecolor="white")
      title("One step prediction error",fontsize=30)
      for iii = 1:numLaps
        ePsiMeas = mpcTraj.closedLoopSEY[2:mpcTraj.idx_end[iii],3,iii]
        ePsiPred = mpcTraj.epsiRef[1:mpcTraj.idx_end[iii]-1,iii]
        err = abs(ePsiMeas - ePsiPred)
        sdata = mpcTraj.closedLoopSEY[1:mpcTraj.idx_end[iii]-1,1,iii]
        plot(sdata,err,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
        hold(true)
      end
      xlabel("s",fontsize=30)
      ylabel(L"err = abs(e_\psi - e_\psi^{Ref}) ",fontsize=30)
      legend()
      grid()


      fig2 = figure(34,facecolor="white")
      errors = zeros(numLaps)
      for iii = 1:numLaps
        ePsiMeas = mpcTraj.closedLoopSEY[2:mpcTraj.idx_end[iii],3,iii]
        ePsiPred = mpcTraj.epsiRef[1:mpcTraj.idx_end[iii]-1,iii]
        err = abs(ePsiMeas - ePsiPred)
        errors[iii] = mean(err)
      end
      plot(1:1:numLaps,errors,"-+",linewidth=2.0,color=getColor(CM))
      xlabel("Lap",fontsize=30)
      ylabel(L"mean(err  = abs(e_\psi - e_\psi^{Ref})) ",fontsize=30)
      legend()
      grid()
    end
    function plotInputs(code::AbstractString)

      # load data from JLD files
      log_path_LMPC     = "$(homedir())/simulations/output-LMPC-$(code).jld"
      d_lmpc          = load(log_path_LMPC)

      # determine how many laps were completed in the file
      mpcTraj    = d_lmpc["mpcTraj"]
      numIter    = d_lmpc["numIter"]


      # determine how many laps were completed in the file
      numLaps = length(find(x -> x>1, numIter))
      if numLaps < 1; error("At least one lap has to be completed succesfully!") end

        CM = ColorManager()
        #toggle_mum(CM)
        # plot data  
        fig = figure(777,facecolor="white")
        title("Inputs",fontsize=30)
        subplot(3,1,1)
        for iii = 1:numLaps
          sdata = mpcTraj.closedLoopSEY[1:mpcTraj.idx_end[iii],1,iii]
          adata = mpcTraj.inputHistory[1:mpcTraj.idx_end[iii],1,iii]
          plot(sdata,adata,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
          hold(true)
        end
        grid(true)
        xlabel(L"s",fontsize=30)
        ylabel(L"a",fontsize=30)
        legend()
        resetColorCycle(CM)

        subplot(3,1,2)
        for iii = 1:numLaps
          sdata = mpcTraj.closedLoopSEY[1:mpcTraj.idx_end[iii],1,iii]
          deltadata = mpcTraj.inputHistory[1:mpcTraj.idx_end[iii],2,iii]
          plot(sdata,deltadata,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
          hold(true)
        end
        grid(true)
        xlabel(L"s",fontsize=30)
        ylabel(L"delta",fontsize=30)
        legend()
        resetColorCycle(CM)

        subplot(3,1,3)
        for iii = 1:numLaps
          sdata = mpcTraj.closedLoopSEY[1:mpcTraj.idx_end[iii],1,iii]
          phidata = mpcTraj.inputHistory[1:mpcTraj.idx_end[iii],3,iii]
          plot(sdata,phidata,"-+",label="Lap $(iii)",linewidth=2.0,color=getColor(CM))
          hold(true)
        end
        grid(true)
        xlabel(L"s",fontsize=30)
        ylabel(L"phi",fontsize=30)
        legend()
        resetColorCycle(CM)
      end

      function animateState2(code::AbstractString,iii::Int64,stateNum::Int64=2,plotOverK::Bool=true,offset::Int64=0,speed::Int64=500 )
        if iii < 3; error("Please specify Iteration from LMPC with iterNum > 2") end
          if offset < 0; error("Offset must be a positive integer") end
            if !(2 <= stateNum <= 5); error("Specify a state number between 2 and 5") end

              #log_path_sim = "$(homedir())/simulations/output-SIM-$(code).jld"
              #log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
              log_path_LMPC   = "$(homedir())/simulations/output-LMPC-$(code).jld"

              #d_rec       = load(log_path_record)
              d_lmpc      = load(log_path_LMPC)



              # from lmpc node
              mpcTraj     = d_lmpc["mpcTraj"] #holds closed loop states
              openLoopSEY = d_lmpc["log_mpc_sol_z"] #open loop mpc states
              numIter     = d_lmpc["numIter"] #holds number of steps per iteration
              coeffConst = d_lmpc["coeffConst"]
              mpcParams   = d_lmpc["mpcParams"]
              mpcParams_pF= d_lmpc["mpcParams_pF"] 

              N_PF = mpcParams_pF.N
              N = mpcParams.N

              # Define array of state names
              stateNames = ["s","ey","epsi","v","rho"]    
              pygui(true)
              fig = figure(facecolor="white")
              PyPlot.grid(true)
              PyPlot.ylabel(stateNames[stateNum],fontsize=20)
              PyPlot.title("Trajectory and open loop for $(stateNames[stateNum])",fontsize=20)

              # TODO: Plotting of best and previous trajectories
              # determine indizes of best and previous lap
              prevInd = mpcTraj.selected_Laps[1,1,iii]
              bestInd = mpcTraj.selected_Laps[1,2,iii]

              # determine length of data points for best trajectory and previous trajectory
              steps_prev = mpcTraj.count[prevInd] - 2
              steps_best = mpcTraj.count[bestInd] - 2
              steps_cl = mpcTraj.idx_end[iii]

              # Define lines to be plotted
              state_prev = plot([],[],"g--",linewidth=2.0,label="Prev C.l.")[1] #Closed loop state trajectory of prev Iteration
              state_best = plot([],[],"m--",linewidth=2.0,label="Best C.l.")[1] #Closed loop state trajectory of best Iteration
              state_ss = plot([],[],"r-",linewidth=2.0,label="Used SS")[1] #Safe set (convex combination (lambda!) of best and prev iteration polynomials )

              state_closedloop = plot([],[],"k-",linewidth=2.0,label="Sys C.L.")[1] # actual closed loop state of current iteration
              #state_mismatch = plot([],[],"c--",linewidth=2.0,label="MPC C.L.")[1] # predicted next closed loop state from mpc
              state_openloop = plot([],[],"b-+",linewidth=2.0,markersize=10.0, label="O.L. Pred")[1] #MPC open loop predictions

              state_approxPrev = plot([],[],"g-",linewidth=2.0, label="SS approx Prev")[1]  #polynomial approximation of prev iteration
              state_approxPointsPrev = plot([],[],"g*",markersize=7.0)[1]                  # used points for polynomial approx
              state_approxBest = plot([],[],"m-",linewidth=2.0, label="SS approx Best")[1]  # polynomial approximation of best iteration
              state_approxPointsBest = plot([],[],"m*",markersize=7.0)[1]                  # used points


              # define position of legend
              PyPlot.legend(fontsize=20,bbox_to_anchor=(1.1, 1.1))

              #set axis limits
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


              ylimBtm = [sMin,eyMin,epsiMin,vMin,rhoMin]
              ylimTop = [sMax,eyMax,epsiMax,vMax,rhoMax]
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

              # TODO: Add texts back about lambda and xf violation
              # print text info about lambda and xF violation 
              # (left, nothing) = axes[:get_xlim]()
              # (nothing,top) = axes[:get_ylim]()
              # top = top - (top/10.0)
              # eps_text = text(left,top,"lambda=",fontsize=15)
              # lambda_text = text(left,top-top/(20.0),"xF Violation=",fontsize=15)


              function init_plot()
                # plot background: previous trajectory (SS) and closed loop trajectory
                if plotOverK

                  state_prev[:set_data]([1:1:steps_prev],[mpcTraj.closedLoopSEY[1:steps_prev,stateNum,prevInd]])
                  state_best[:set_data]([1:1:steps_best],[mpcTraj.closedLoopSEY[1:steps_best,stateNum,bestInd]])
                  state_ss[:set_data]([(N+1):1:steps_cl+N],[mpcTraj.xfStates[1:steps_cl,stateNum,iii]])
                  state_closedloop[:set_data]([1:1:steps_cl],[mpcTraj.closedLoopSEY[1:steps_cl,stateNum,iii]])
                  #state_mismatch[:set_data]([2:1:steps_cl],[Iterations[iii].openLoopSEY[stateNum,2,1:steps_cl-1]]) #FIXME Also check here if ytrue time index used or it has to be t-1
                else
                  state_prev[:set_data]([mpcTraj.closedLoopSEY[1:steps_prev,1,prevInd]],[mpcTraj.closedLoopSEY[1:steps_prev,stateNum,prevInd]])
                  state_best[:set_data]([mpcTraj.closedLoopSEY[1:steps_best,1,bestInd]],[mpcTraj.closedLoopSEY[1:steps_best,stateNum,bestInd]])
                  state_ss[:set_data]([mpcTraj.xfStates[1:steps_cl,1,iii]],[mpcTraj.xfStates[1:steps_cl,stateNum,iii]]) 
                  state_closedloop[:set_data]([mpcTraj.closedLoopSEY[1:steps_cl,1,iii]],[mpcTraj.closedLoopSEY[1:steps_cl,stateNum,iii]])
                  #state_mismatch[:set_data]([Iterations[iii].openLoopSEY[1,2,:]],[Iterations[iii].openLoopSEY[stateNum,2,:]])
                end

                return (state_prev,state_best,state_ss,state_closedloop,nothing)
              end

              function animate(frameNum)
                frameNum += 1
                k = frameNum + offset


                ss_ind = [prevInd,bestInd]
                approxLines = [state_approxPointsPrev state_approxPointsBest; state_approxPrev state_approxBest]

                if plotOverK
                  # Update open loop trajectory data
                  state_openloop[:set_data]([k:1:k+N],[openLoopSEY[k,1:N+1,stateNum,iii]])
                  for i=1:2
                    ind = ss_ind[i]
                    xfRangeMin = Int64(mpcTraj.xfRange[k,1,i,iii])
                    xfRangeMax = Int64(mpcTraj.xfRange[k,2,i,iii])
                    sRange = mpcTraj.closedLoopSEY[xfRangeMin:xfRangeMax,1,ind]
                    stateRange = mpcTraj.closedLoopSEY[xfRangeMin:xfRangeMax,stateNum,ind]

                    # plot approximation points used for interpolation over k
                    approxLines[1,i][:set_data]([xfRangeMin:1:xfRangeMax],[stateRange])


                    # get polynomail coefficients
                    ss = coeffConst[:,i,stateNum-1,k,iii]
                    # determine used polynomial order
                    poly_order = size(ss,1) - 1 
                    # evaluate polynomial in sRange
                    statePoly = sum(sRange.^(poly_order+1-p) * ss[p] for p=1:poly_order + 1 )
                    approxLines[2,i][:set_data]([xfRangeMin:1:xfRangeMax],[statePoly])
                  end      
                else
                  # Update open loop trajectory data
                  state_openloop[:set_data]([openLoopSEY[k,1:N+1,1,iii]],[openLoopSEY[k,1:N+1,stateNum,iii]])

                  for i=1:2
                    ind = ss_ind[i]
                    xfRangeMin = Int64(mpcTraj.xfRange[k,1,i,iii])
                    xfRangeMax = Int64(mpcTraj.xfRange[k,2,i,iii])
                    sRange = mpcTraj.closedLoopSEY[xfRangeMin:xfRangeMax,1,ind]
                    stateRange = mpcTraj.closedLoopSEY[xfRangeMin:xfRangeMax,stateNum,ind]

                    # plot approximation points used for interpolation over k
                    approxLines[1,i][:set_data]([sRange],[stateRange])

                    # get polynomail coefficients
                    ss = coeffConst[:,i,stateNum-1,k,iii]
                    # evaluate polynomial in sRange
                    poly_order = size(ss,1) - 1 
                    statePoly = sum(sRange.^(poly_order+1-p) * ss[p] for p=1:poly_order + 1 )
                    approxLines[2,i][:set_data]([sRange],[statePoly])
                  end
                end

                # # update text about softconstrain violation and color
                # lambda_text[:set_text]("lambda = $(Iterations[iii].lambdas[k])")
                # eps_text[:set_text]("xF Violation = $(Iterations[iii].epsvec[stateNum-1,k])")
                # if Iterations[iii].epsvec[stateNum-1,k] > 0.01
                #   eps_text[:set_color]("r")
                # else
                #   eps_text[:set_color]("g")
                # end

                # return (state_openloop,state_approxPrev,state_approxPointsPrev,state_approxBest, state_approxPointsBest,eps_text,lambda_text,nothing)
                return (state_openloop,nothing)
              end
              NumFrames = length(mpcTraj.closedLoopSEY[:,1,iii]) - offset
              anim = animation.FuncAnimation(fig, animate, init_func=init_plot, frames=NumFrames, interval=speed,repeat=false)
            end





            function showCost(code::AbstractString)
              log_path_LMPC   = "$(homedir())/simulations/output-LMPC-$(code).jld"
              log_path_record = "$(homedir())/simulations/output-record-$(code).jld"
              d_rec       = load(log_path_record)
              d_lmpc      = load(log_path_LMPC)

              t           = d_lmpc["t"]
              state       = d_lmpc["state"]
              cost        = d_lmpc["cost"]


              t0 = t[1]

              figure()
              ax1=subplot(211)
              title("MPC states and cost")
              plot(t-t0,state)
              legend(["s","ey","epsi","v","rho","epsiRef","filter"])
              grid(1)
              subplot(212,sharex = ax1)
              plot(t-t0,cost)
              grid(1)
              legend(["costZ","costZTerm (term. constr.)","constZTerm (term. cost)","derivCost","controlCost","laneCost","modelErrorCost"])
            end



            function create_track(w)
              x = [0.0]           # starting point
              y = [0.0]
              x_l = [0.0]           # starting point
              y_l = [w]
              x_r = [0.0]           # starting point
              y_r = [-w]
              ds = 0.06

              theta = [0.0]

              # SOPHISTICATED TRACK
              # add_curve(theta,30,0.0)
              # add_curve(theta,60,-2*pi/3)
              # add_curve(theta,90,pi)
              # add_curve(theta,80,-5*pi/6)
              # add_curve(theta,10,0.0)
              # add_curve(theta,50,-pi/2)
              # add_curve(theta,50,0.0)
              # add_curve(theta,40,-pi/4)
              # add_curve(theta,30,pi/4)
              # add_curve(theta,20,0.0)
              # add_curve(theta,50,-pi/2)
              # add_curve(theta,25,0.0)
              # add_curve(theta,50,-pi/2)
              # add_curve(theta,28,0.0)

              # # SIMPLE track
              # add_curve(theta,50,0)
              # add_curve(theta,100,-pi)
              # add_curve(theta,100,0)
              # add_curve(theta,100,-pi)
              # add_curve(theta,49,0)

              # GOGGLE TRACK
              # add_curve(theta,30,0)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,20,-pi/6)
              # add_curve(theta,30,pi/3)
              # add_curve(theta,20,-pi/6)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,35,0)

              # SIMPLE GOGGLE TRACK
              # add_curve(theta,30,0)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,10,0)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,20,pi/10)
              # add_curve(theta,30,-pi/5)
              # add_curve(theta,20,pi/10)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,10,0)
              # add_curve(theta,40,-pi/2)
              # add_curve(theta,35,0)

              # MICHAEL'S TRACK
              add_curve(theta,25,0)
              add_curve(theta,40,-pi/2)
              add_curve(theta,45,0)
              add_curve(theta,40,-pi/2)
              add_curve(theta,10,0)
              add_curve(theta,40,-pi/2)
              add_curve(theta,5,0)
              add_curve(theta,40,+pi/2)
              add_curve(theta,6,0)
              add_curve(theta,30,-pi/2)
              add_curve(theta,5,0)
              add_curve(theta,30,-pi/2)
              add_curve(theta,38,0)

              #  # SHORT SIMPLE track
              # add_curve(theta,10,0)
              # add_curve(theta,80,-pi)
              # add_curve(theta,20,0)
              # add_curve(theta,80,-pi)
              # add_curve(theta,9,0)

              for i=1:length(theta)
                push!(x, x[end] + cos(theta[i])*ds)
                push!(y, y[end] + sin(theta[i])*ds)
                push!(x_l, x[end-1] + cos(theta[i]+pi/2)*w)
                push!(y_l, y[end-1] + sin(theta[i]+pi/2)*w)
                push!(x_r, x[end-1] + cos(theta[i]-pi/2)*w)
                push!(y_r, y[end-1] + sin(theta[i]-pi/2)*w)
              end
              return(x, y, x_l, y_l, x_r, y_r)
            end

            function add_curve(theta::Array{Float64}, length::Int64, angle)
              d_theta = 0
              curve = 2*sum(1:length/2)+length/2
              for i=0:length-1
                if i < length/2+1
                  d_theta = d_theta + angle / curve
                else
                  d_theta = d_theta - angle / curve
                end
                push!(theta, theta[end] + d_theta)
              end
            end

          end #MODULE
          # add_curve(theta,30,0)
    # add_curve(theta,20,-pi/6)