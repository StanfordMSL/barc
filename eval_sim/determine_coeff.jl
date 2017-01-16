    
using JLD, ProfileView
    code = "59b7"
    log_path_LMPC   = "$(homedir())/simulations/coefficients_charged/output-LMPC-$(code).jld"
    d_lmpc      = load(log_path_LMPC)

    c_Vx        = d_lmpc["c_Vx"]
    c_Vy        = d_lmpc["c_Vy"]
    c_Psi       = d_lmpc["c_Psi"]


    relevant_data_start = findfirst(f->!isnan(f),c_Vx) +10
   
    last = 700
    c_vx_1_rel = c_Vx[relevant_data_start:end,1]
    c_vx_2_rel = c_Vx[relevant_data_start:end,2]
    c_vx_3_rel = c_Vx[relevant_data_start:end,3]

    c_vy_1_rel = c_Vy[relevant_data_start:end,1]
    c_vy_2_rel = c_Vy[relevant_data_start:end,2]
    c_vy_3_rel = c_Vy[relevant_data_start:end,3]
    c_vy_4_rel = c_Vy[relevant_data_start:end,4]

    c_psi_1_rel = c_Psi[relevant_data_start:end,1]
    c_psi_2_rel = c_Psi[relevant_data_start:end,2]
    c_psi_3_rel = c_Psi[relevant_data_start:end,3]

    println("------------------ Calculate averages of coefficent data: ------------------\n")
    println("mean(C_Vx[1]): $(mean(c_vx_1_rel))\n")
    println("mean(C_Vx[2]): $(mean(c_vx_2_rel))\n")
    println("mean(C_Vx[3]): $(mean(c_vx_3_rel))\n \n")

    println("mean(C_Vy[1]): $(mean(c_vy_1_rel))\n")
    println("mean(C_Vy[2]): $(mean(c_vy_2_rel))\n")
    println("mean(C_Vy[3]): $(mean(c_vy_3_rel))\n")
    println("mean(C_Vy[4]): $(mean(c_vy_4_rel))\n \n")

    println("mean(C_Psi[1]): $(mean(c_psi_1_rel))\n")
    println("mean(C_Psi[2]): $(mean(c_psi_2_rel))\n")
    println("mean(C_Psi[3]): $(mean(c_psi_3_rel))\n")