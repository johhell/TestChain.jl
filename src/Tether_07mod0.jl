#=
   This code is a modified version of an algorithm
   forming part of the software package
   Tethers https://github.com/ufechner7/Tethers.jl
   authored by  Uwe Fechner and Bart van de Lint.
   The original software code is subject to the MIT License.
=#


module T07mod0

# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots

@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 50                                      # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 5                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 10                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
function calc_initial_state(se)
    POS0 = zeros(3, se.segments+1)
    VEL0 = zeros(3, se.segments+1)
    for i in 1:se.segments+1
        l0 = -(i-1)*se.l0/se.segments
        v0 = -(i-1)*se.v_ro/se.segments
        POS0[:, i] .= [sin(se.α0) * l0, 0, cos(se.α0) * l0]
        VEL0[:, i] .= [sin(se.α0) * v0, 0, cos(se.α0) * v0]
    end
    POS0, VEL0
end

function model(se)
    POS0, VEL0 = calc_initial_state(se)
    mass_per_meter = se.rho_tether * π * (se.d_tether/2000.0)^2
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin 
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[1:se.segments]
        rel_vel(t)[1:3, 1:se.segments]
        spring_vel(t)[1:se.segments]
        c_spr(t)[1:se.segments]
        spring_force(t)[1:3, 1:se.segments]
        v_apparent(t)[1:3, 1:se.segments]
        v_app_perp(t)[1:3, 1:se.segments]
        norm_v_app(t)[1:se.segments]
        half_drag_force(t)[1:3, 1:se.segments]
        total_force(t)[1:3, 1:se.segments+1]
    end
    # basic differential equations
    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    # loop over all segments to calculate the spring and drag forces
    for i in 1:se.segments
        eqs = [segment[:, i]      ~ pos[:, i+1] - pos[:, i],
               len[i]             ~ norm(segment[:, i]),
               unit_vector[:, i]  ~ -segment[:, i]/len[i],
               rel_vel[:, i]      ~ vel[:, i+1] - vel[:, i],
               spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
               c_spr[i]           ~ c_spring / (1+rel_compression_stiffness) 
                                     * (rel_compression_stiffness+(len[i] > l_spring)),
               spring_force[:, i] ~ (c_spr[i] * (len[i] - l_spring) 
                                     + damping * spring_vel[i]) * unit_vector[:, i],
               v_apparent[:, i]   ~ se.v_wind_tether .- (vel[:, i] + vel[:, i+1])/2,
               v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
               norm_v_app[i]      ~ norm(v_app_perp[:, i]),
               half_drag_force[:, i] ~ 0.25 * se.rho * se.cd_tether * norm_v_app[i] * (len[i]*se.d_tether/1000.0) * v_app_perp[:, i]
            ]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # loop over all tether particles to apply the forces and calculate the accelerations
    for i in 1:(se.segments+1)
        eqs = []
        if i == se.segments+1
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] + half_drag_force[:, i-1])
#             push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / (0.5 * m_tether_particle))  # ???
            push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / ( m_tether_particle))
        elseif i == 1
            push!(eqs, total_force[:, i] ~ spring_force[:, i] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ zeros(3))
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1] - spring_force[:, i] 
                                           + half_drag_force[:, i-1] + half_drag_force[:, i])
            push!(eqs, acc[:, i]         ~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # scalar equations
    eqs = [l_spring          ~ (se.l0 + se.v_ro*t) / se.segments,
           c_spring          ~ se.c_spring / l_spring,
           m_tether_particle ~ mass_per_meter * l_spring,
           damping           ~ se.damping  / l_spring]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))  
        
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel, l_spring, c_spr
end

function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end


se = Settings3(; segments=5,v_ro=0.0, v_wind_tether=[0.00,0.0,0.01], rho=0.0)
set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
simple_sys, pos, vel, len, c_spr = model(se)
sol, elapsed_time = simulate(se, simple_sys)




T = sol.t


using HDF5

h5open("T07mod.hdf5", "w") do fid
    fid["zeit"] = T
    fid["model"] = "T07mod-MTK0"
    fid["segments"] = se.segments
    gruppe= create_group(fid, "positions")
    for i = 1:se.segments
        xx = sol(T, idxs=pos[1,i+1]).u
        yy = sol(T, idxs=pos[2,i+1]).u
        zz = sol(T, idxs=pos[3,i+1]).u
        gruppe["pos$(i)"] = hcat(xx,yy,zz)
    end
end



nothing

end
