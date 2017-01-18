type MpcModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    coeffTermConst::Array{JuMP.NonlinearParameter,3}
    coeffTermCost::Array{JuMP.NonlinearParameter,2}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    ParInt::JuMP.Variable

    laneCost::JuMP.NonlinearExpression
    constZTerm::JuMP.NonlinearExpression
    costZTerm::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression

    uPrev::Array{JuMP.NonlinearParameter,2}

    function MpcModel(mpcParams::MpcParams,mpcCoeff::MpcCoeff,modelParams::ModelParams,trackCoeff::TrackCoeff)
        m = new()
        dt              = modelParams.dt
        L_a             = modelParams.l_A
        L_b             = modelParams.l_B
        barcWidth       = modelParams.width
        c0              = modelParams.c0
        c_Vx            = modelParams.c_Vx::Array{Float64,1}
        c_Vy            = modelParams.c_Vy::Array{Float64,1}
        c_Psi            = modelParams.c_Psi::Array{Float64,1}


        N               = mpcParams.N
        Q               = mpcParams.Q
        Q_term          = mpcParams.Q_term
        R               = mpcParams.R
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        ey_max          = trackCoeff.width/2 - barcWidth/2 #take into account the vehicle width

        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        Q_term_cost     = mpcParams.Q_term_cost::Float64
        delay_df        = mpcParams.delay_df
        delay_a         = mpcParams.delay_a
        u_lb            = mpcParams.u_lb::Array{Float64,1}
        u_ub            = mpcParams.u_ub::Array{Float64,1}
        z_lb            = mpcParams.z_lb::Array{Float64,1}
        z_ub            = mpcParams.z_ub::Array{Float64,1}
        du_lb           = mpcParams.du_lb::Array{Float64,1}  #bounds on hard input derivative constraints
        du_ub           = mpcParams.du_ub::Array{Float64,1}

        acc_f           = 1.0


        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

        mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=0.08))

        @variable( mdl, z_Ol[1:(N+1),1:7])
        @variable( mdl, u_Ol[1:N,1:2])
        @variable( mdl, 0 <= ParInt <= 1)
        @variable( mdl, eps[1:N+1] >= 0) # eps for soft lane constraints

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb[i])
                setupperbound(u_Ol[j,i], u_ub[i])
            end
        end
        for i=1:7
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb[i])
                setupperbound(z_Ol[j,i], z_ub[i])
            end
        end

        @NLparameter(mdl, z0[i=1:7] == 0)
        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)
        @NLparameter(mdl, coeffTermConst[i=1:order+1,j=1:2,k=1:5] == 0)
        @NLparameter(mdl, coeffTermCost[i=1:order+1,j=1:2] == 0)
        @NLparameter(mdl, uPrev[1:10,1:2] == 0)

        # Conditions for first solve:
        setvalue(z0[1],1)

        @NLconstraint(mdl, [i=1:7], z_Ol[1,i] == z0[i])
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,5] <= ey_max + eps[i])
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,5] >= -ey_max - eps[i])

        # Expressions to formulate the system dynamics
        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,6]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])
        @NLexpression(mdl, dsdt[i = 1:N], (z_Ol[i,1]*cos(z_Ol[i,4]) - z_Ol[i,2]*sin(z_Ol[i,4]))/(1-z_Ol[i,5]*c[i]))
        
        println("Initializing model...")

        # System dynamics
        for i=1:N
            if i<=delay_df
                @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*uPrev[delay_df+1-i,2]) # yDot
                @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] + c_Psi[3]*uPrev[delay_df+1-i,2])                            # psiDot
            else
                @NLconstraint(mdl, z_Ol[i+1,2]  == z_Ol[i,2] + c_Vy[1]*z_Ol[i,2]/z_Ol[i,1] + c_Vy[2]*z_Ol[i,1]*z_Ol[i,3] + c_Vy[3]*z_Ol[i,3]/z_Ol[i,1] + c_Vy[4]*u_Ol[i-delay_df,2]) # yDot
                @NLconstraint(mdl, z_Ol[i+1,3]  == z_Ol[i,3] + c_Psi[1]*z_Ol[i,3]/z_Ol[i,1] + c_Psi[2]*z_Ol[i,2]/z_Ol[i,1] + c_Psi[3]*u_Ol[i-delay_df,2])                            # psiDot
            end
            if i<=delay_a
                @NLconstraint(mdl, z_Ol[i+1,7]  == z_Ol[i,7] + dt*(uPrev[delay_a+1-i,1]-z_Ol[i,7])*acc_f)
            else
                @NLconstraint(mdl, z_Ol[i+1,7]  == z_Ol[i,7] + dt*(u_Ol[i-delay_a,1]-z_Ol[i,7])*acc_f)
            end

            @NLconstraint(mdl, z_Ol[i+1,1]  == z_Ol[i,1] + c_Vx[1]*z_Ol[i,2]*z_Ol[i,3] + c_Vx[2]*z_Ol[i,1] + c_Vx[3]*z_Ol[i,7]) 
            @NLconstraint(mdl, z_Ol[i+1,4]  == z_Ol[i,4] + dt*(z_Ol[i,3]-dsdt[i]*c[i]))                                                                                 # ePsi
            @NLconstraint(mdl, z_Ol[i+1,5]  == z_Ol[i,5] + dt*(z_Ol[i,1]*sin(z_Ol[i,4])+z_Ol[i,2]*cos(z_Ol[i,4])))                                                      # eY
            @NLconstraint(mdl, z_Ol[i+1,6]  == z_Ol[i,6] + dt*dsdt[i]  )                                                                                                # s
        end

        # Hard input derivative constraints
        @NLconstraint(mdl, u_Ol[1,2]-uPrev[1,2] <= du_ub[2])
        @NLconstraint(mdl, u_Ol[1,2]-uPrev[1,2] >= du_lb[2])

        #@NLconstraint(mdl, u_Ol[1,1]-uPrev[1,1] <= du_ub[1])
        #@NLconstraint(mdl, u_Ol[1,1]-uPrev[1,1] >= du_lb[1])
        #for j=1:2
            for i=1:N-1 # Constraints on u:
                @NLconstraint(mdl, u_Ol[i+1,2]-u_Ol[i,2] <= du_ub[2])
                @NLconstraint(mdl, u_Ol[i+1,2]-u_Ol[i,2] >= du_lb[2])
            end
        #end

        # Cost functions

        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:6} +
                                          QderivU[1]*((uPrev[1,1]-u_Ol[1,1])^2+sum{(u_Ol[i,1]-u_Ol[i+1,1])^2,i=1:N-delay_a-1})+
                                          QderivU[2]*((uPrev[1,2]-u_Ol[1,2])^2+sum{(u_Ol[i,2]-u_Ol[i+1,2])^2,i=1:N-delay_df-1}))

        # Lane cost
        # ---------------------------------
        #@NLexpression(mdl, laneCost, sum{100000*eps[i]+1000*eps[i]^2,i=1:2})
        @NLexpression(mdl, laneCost, sum{10*eps[i]+100*eps[i]^2,i=2:N+1})
        
        # Lane cost
        # ---------------------------------
        #@NLexpression(mdl, laneCost, 100*sum{z_Ol[i,5]^2*((0.5+0.5*tanh(10*(z_Ol[i,5]-ey_max))) + (0.5-0.5*tanh(10*(z_Ol[i,5]+ey_max)))),i=1:N+1})

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, R[1]*sum{(u_Ol[i,1])^2,i=1:N-delay_a}+
                                        R[2]*sum{(u_Ol[i,2])^2,i=1:N-delay_df})

        # Terminal constraints (soft), starting from 2nd lap
        # ---------------------------------
        @NLexpression(mdl, constZTerm, sum{Q_term[j]*(ParInt*(sum{coeffTermConst[i,1,j]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermConst[order+1,1,j])+
                                            (1-ParInt)*(sum{coeffTermConst[i,2,j]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermConst[order+1,2,j])-z_Ol[N+1,j])^2,j=1:5})
        
        # Terminal cost
        # ---------------------------------
        # The value of this cost determines how fast the algorithm learns. The higher this cost, the faster the control tries to reach the finish line.
        @NLexpression(mdl, costZTerm, 1/5*(Q_term_cost*(ParInt*(sum{coeffTermCost[i,1]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermCost[order+1,1])+
                                      (1-ParInt)*(sum{coeffTermCost[i,2]*z_Ol[N+1,6]^(order+1-i),i=1:order}+coeffTermCost[order+1,2]))))
        

        @NLexpression(mdl, costZ, 0*QderivU[1]*sum{z_Ol[i,7],i=1:N})
        # Solve model once
        @NLobjective(mdl, Min, derivCost + constZTerm + costZTerm + laneCost + N*Q_term_cost)
        #@NLobjective(mdl, Min, derivCost + costZ)
        sol_stat=solve(mdl)
        println("Finished solve 1: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2: $sol_stat")
        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.ParInt = ParInt
        m.uPrev = uPrev

        m.coeffTermCost = coeffTermCost
        m.coeffTermConst = coeffTermConst

        m.derivCost = derivCost
        m.controlCost = controlCost
        m.laneCost = laneCost
        m.constZTerm = constZTerm
        m.costZTerm  = costZTerm
        m.costZ  = costZ
        return m
    end
end

type MpcModel_pF
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}

    derivCost::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    uPrev::Array{JuMP.NonlinearParameter,2}

    function MpcModel_pF(mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff)
        println("Starting creation of pf model")
        m = new()
        dt          = modelParams.dt
        L_a         = modelParams.l_A
        L_b         = modelParams.l_B
        c0          = modelParams.c0
  

        N           = mpcParams.N
        Q           = mpcParams.Q
        R           = mpcParams.R
        QderivZ     = mpcParams.QderivZ::Array{Float64,1}
        QderivU     = mpcParams.QderivU::Array{Float64,1}
        delay_df    = mpcParams.delay_df::Int64
        delay_a     = mpcParams.delay_a::Int64
        u_lb        = mpcParams.u_lb::Array{Float64,1}
        u_ub        = mpcParams.u_ub::Array{Float64,1}
        z_lb        = mpcParams.z_lb::Array{Float64,1}
        z_ub        = mpcParams.z_ub::Array{Float64,1}

        v_ref       = mpcParams.vPathFollowing

        acc_f       = 1.0

        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

        # Create function-specific parameters
        z_Ref::Array{Float64,2}
        z_Ref       = cat(2,zeros(N+1,3),v_ref*ones(N+1,1))       # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref       = zeros(N,2)

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=0.07))

        # Create variables (these are going to be optimized)
        @variable( mdl, z_Ol[1:(N+1),1:5], start = 0)          # z = s, ey, epsi, v
        @variable( mdl, u_Ol[1:N,1:2], start = 0)

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb[i])
                setupperbound(u_Ol[j,i], u_ub[i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb[i])
                setupperbound(z_Ol[j,i], z_ub[i])
            end
        end

        @NLparameter(mdl, z0[i=1:5] == 0)
        @NLparameter(mdl, uPrev[1:10,1:2] == 0)

        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == 0)

        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1])

        # System dynamics
        setvalue(z0[4],v_ref)
        @NLconstraint(mdl, [i=1:5], z_Ol[1,i] == z0[i])         # initial condition
        for i=1:N
            if i<=delay_df
                @NLexpression(mdl, bta[i],  atan( L_a / (L_a + L_b) * tan( uPrev[delay_df+1-i,2] ) ) )
            else
                @NLexpression(mdl, bta[i],  atan( L_a / (L_a + L_b) * tan( u_Ol[i-delay_df,2] ) ) )
            end
            if i<=delay_a
                #@NLconstraint(mdl, z_Ol[i+1,4] == z_Ol[i,4] + dt*(uPrev[delay_a+1-i,1] - 0.5*z_Ol[i,4]))  # v
                @NLconstraint(mdl, z_Ol[i+1,5] == z_Ol[i,5] + dt*(uPrev[delay_a+1-i,1] - z_Ol[i,5])*acc_f)  # v
            else
                #@NLconstraint(mdl, z_Ol[i+1,4] == z_Ol[i,4] + dt*(u_Ol[i-delay_a,1] - 0.5*z_Ol[i,4]))     # v
                @NLconstraint(mdl, z_Ol[i+1,5] == z_Ol[i,5] + dt*(u_Ol[i-delay_a,1] - z_Ol[i,5])*acc_f)     # v
            end

            @NLexpression(mdl, dsdt[i], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))
            @NLconstraint(mdl, z_Ol[i+1,1] == z_Ol[i,1] + dt*dsdt[i]  )                                                # s
            @NLconstraint(mdl, z_Ol[i+1,2] == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                        # ey
            @NLconstraint(mdl, z_Ol[i+1,3] == z_Ol[i,3] + dt*(z_Ol[i,4]/L_a*sin(bta[i])-dsdt[i]*c[i])  )               # epsi
            @NLconstraint(mdl, z_Ol[i+1,4] == z_Ol[i,4] + dt*(z_Ol[i,5] - 0.5*z_Ol[i,4]))  # v
        end

        # Cost definitions
        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:4} +
                                            QderivU[1]*((uPrev[1,1]-u_Ol[1,1])^2+sum{(u_Ol[i,1]-u_Ol[i+1,1])^2,i=1:N-delay_a-1})+
                                            QderivU[2]*((uPrev[1,2]-u_Ol[1,2])^2+sum{(u_Ol[i,2]-u_Ol[i+1,2])^2,i=1:N-delay_df-1}))

        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*R[1]*sum{(u_Ol[i,1])^2,i=1:N-delay_a}+
                                        0.5*R[2]*sum{(u_Ol[i,2])^2,i=1:N-delay_df})

        # State cost
        # ---------------------------------
        @NLexpression(mdl, costZ, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=2:N+1},i=1:4})    # Follow trajectory

        # Objective function
        @NLobjective(mdl, Min, costZ + derivCost + controlCost)

        # create first artificial solution (for warm start)
        for i=1:N+1
            setvalue(z_Ol[i,:],[(i-1)*dt*v_ref 0 0 v_ref 0])
        end
        for i=1:N
            setvalue(u_Ol[i,:],[0.15 0])
        end
        # First solve
        sol_stat=solve(mdl)
        println("Finished solve 1: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2: $sol_stat")
        
        m.mdl = mdl
        m.z0 = z0
        m.coeff = coeff
        m.z_Ol = z_Ol
        m.u_Ol = u_Ol
        m.uPrev = uPrev
        m.derivCost = derivCost
        m.costZ = costZ
        m.controlCost = controlCost
        return m
    end
end

# why does the solution change when I start the simulation twice, totally independently? Is not everything initialized equally?
# why do I get a failed restoration phase even though there are *almost* no constraints and the solution should be super clear?

# does it actually make sense to use the first (initial) state as a @variable and combine it with a constraint to fix it to z0 ?
# or better use only a parameter for z0? -> Found out that solutions are different if I set constraints for first state or not!?
