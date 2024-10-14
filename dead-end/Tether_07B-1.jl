# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (no spring forces
# for l < l_0), n tether segments and reel-in and reel-out. 


using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots


include("./Tparameters.jl")


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
#     @parameters rel_compression_stiffness = se.rel_compression_stiffness


    @variables begin
        pos(t)[1:3, 1:se.segments+1]  = POS0
        vel(t)[1:3, 1:se.segments+1]  = VEL0
        acc(t)[1:3, 1:se.segments+1]
        segment(t)[1:3, 1:se.segments]
        unit_vector(t)[1:3, 1:se.segments]
        length(t), c_spring(t), damping(t), m_tether_particle(t)
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


    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    for i in se.segments:-1:1
        eqs = [ segment[:, i] ~ pos[:, i+1] - pos[:, i],
            len[i] ~ norm(segment[:, i]),
            unit_vector[:, i] ~ -segment[:, i]/len[i],
            rel_vel[:, i] ~ vel[:, i+1] - vel[:, i],
            spring_vel[i] ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
            c_spr[i] ~ c_spring * (len[i] > length/se.segments),
            spring_force[:, i] ~ (c_spr[i] * (len[i] - (length/se.segments)) + damping * spring_vel[i]) * unit_vector[:, i],
            ]
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end


    for i in 1:(se.segments+1)
        eqs = []
        if i == 1   #fist node
            push!(eqs, total_force[:, 1] ~ spring_force[:, 1]) # forces are applied, but fixed position
            push!(eqs, acc[:, 1] ~ zeros(3))    #FIXED position
        elseif i == (se.segments+1) #letzter Knoten
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1])
            push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / m_tether_particle)
        else
            push!(eqs, total_force[:, i] ~ spring_force[:, i-1]- spring_force[:, i])
            push!(eqs, acc[:, i] ~ se.g_earth + total_force[:, i] / m_tether_particle)
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end

    eqs = [acc[:, 1] ~ zeros(3),
            length ~ se.l0 + se.v_ro*t,
            c_spring ~ se.c_spring / (length/se.segments),
            m_tether_particle ~ mass_per_meter * (length/se.segments),
            damping  ~ se.damping  / (length/se.segments),  ]
    eqs2 = vcat(eqs2, reduce(vcat, eqs))

    Syssys = reduce(vcat, Symbolics.scalarize.(eqs2))
    @named sys = ODESystem(Syssys, t)
    simple_sys = structural_simplify(sys)
    simple_sys, pos, vel, c_spr
end



function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    @time sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)
    sol
end

function plot2d(se, sol, pos, reltime, line, sc, txt, j)
    index = Int64(round(reltime*50+1))
    x, z = Float64[], Float64[]
    for particle in 1:se.segments+1
        push!(x, (sol(sol.t, idxs=pos[1, particle]))[index])
        push!(z, (sol(sol.t, idxs=pos[3, particle]))[index])
    end
    z_max = maximum(z)
    if isnothing(line)
        line, = plot(x,z; linewidth="1")
        sc  = scatter(x, z; s=15, color="red") 
        txt = annotate("t=$(round(reltime,digits=1)) s",  
                        xy=(se.l0/4.2, z_max-7), fontsize = 12)
    else
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        gcf().canvas.draw()
    end
    if se.save
        PyPlot.savefig("video/"*"img-"*lpad(j,4,"0"))
    end
    line, sc, txt
end




se = Settings(; v_ro = 0.0)   # reel-out speed = 0

simple_sys, pos, vel, c_spr = model(se)

sol = simulate(se, simple_sys)

pos = sol(sol.t, idxs=pos)

nothing

#=
julia> pos[end]
3×3 Matrix{Float64}:
  1.77981e-15   -0.943194   -0.645232
  0.0            0.0         0.0
 -2.66907e-17  -24.9824    -49.9807

B5modia
 -0.0581835  0.0  -49.9869
  -0.204829   0.0  -49.986
  -0.351559   0.0  -49.9846
  -0.498363   0.0  -49.9829
  -0.645229   0.0  -49.9807

=#
