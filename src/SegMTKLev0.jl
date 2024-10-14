
module CH8seg

using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit
using DifferentialEquations
using LinearAlgebra


include("Tparameters.jl")


function Segment(param::Settings, idx::Int64; name)

    @variables begin
        pos(t)[1:3]
        posM(t)[1:3] = param.positionStart[idx]
        velM(t)[1:3] = [0.0,0.0,0.0]
        accM(t)[1:3]
        segment(t)[1:3]
        segUnit(t)[1:3]
        vRel(t)[1:3]
        vRel_unit(t)[1:3]
        forceFriction(t)[1:3]
        forceSpring(t)[1:3]
        forceExternal(t)[1:3]
        segmentABS(t)
        vRelABS(t)
        v_seg(t)
        forceSpringABS(t)
        lgwind(t)
    end

    eqs = [
        D(posM) ~ velM,
        D(velM) ~ accM,
        segment ~ posM - pos,
        segUnit ~ segment / segmentABS,
        vRel ~ velM - param.wind_speed,
        vRelABS ~ max(0.01,norm(vRel)),
        vRel_unit ~ vRel / vRelABS,
        lgwind ~ norm(cross(vRel_unit,segment)),
        forceFriction ~ lgwind * vRel_unit * param.frictionCoeff * vRelABS^2,
        segmentABS ~ norm(segment),
        v_seg ~ D(segmentABS),
        forceSpringABS ~ -(segmentABS-param.len_per_segment)*param.c_spring0 - v_seg*param.damping,
        forceSpring ~ segUnit * forceSpringABS,
        (accM - param.g_earth)* param.mass_per_seg ~  forceSpring + forceExternal - forceFriction,
        ]
    eqs2 = reduce(vcat, eqs)
    Syssys = reduce(vcat, Symbolics.scalarize.(eqs2))
    ODESystem(Syssys, t; name = name)
end



listofSegments = []

function AddNewSegment(idx::Int, segList::Vector, ss::Settings)
    s = Symbol("S", idx)
    @eval ( @named ($s) = Segment($ss,$idx))
    @eval ( push!($segList,($s)))
end


function model(se::Settings)

    # creating segments
    for i = 1:se.segments
        AddNewSegment(i, listofSegments, se)
    end

    # connecting segments
    #TODO   Future: could be replaced by @connnector???
    function SegmentConnections(listSegs::Vector)
        eqs = []
        eqs = vcat(eqs, listSegs[1].pos ~ [0.0,0.0, 0.0])  # 1st segment with fixed position
        for i = 1:se.segments-1  # connecting segments
            eqs = vcat(eqs, listSegs[i+1].pos ~ listSegs[i].posM)   # position
            eqs = vcat(eqs, listSegs[i].forceExternal ~ - listSegs[i+1].forceSpring) # force
        end
        eqs = vcat(eqs, listSegs[end].forceExternal ~ se.endForce)     # force @last segment
    end

    connectEqu = SegmentConnections(listofSegments)
    @named sys = ODESystem(connectEqu, t, systems=listofSegments)
    return simpsys = structural_simplify(sys)
end



function simulate(se, simple_sys)
    dt = 0.02
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    prob = ODEProblem(simple_sys, nothing, tspan)
    sol = solve(prob, Rodas5(), dt=dt, abstol=tol, reltol=tol, saveat=ts)
end
#

se = Settings(; v_ro = 0.0, frictionCoeff = 0.0, segments=5)

simple_sys = model(se)

sol = simulate(se, simple_sys)



# position of the last particle


T = sol.t

using HDF5

h5open("SegMTKL0.hdf5", "w") do fid
    fid["zeit"] = T
    fid["model"] = "Seg-MTK0"
    fid["segments"] = se.segments
    gruppe= create_group(fid, "positions")
    for i = 1:se.segments
        E = listofSegments[i]
        xx = sol(T, idxs=E.posM[1]).u
        yy = sol(T, idxs=E.posM[2]).u
        zz = sol(T, idxs=E.posM[3]).u
        gruppe["pos$(i)"] = hcat(xx,yy,zz)
    end
end



using Plots
lastSeg = listofSegments[end]
xx = sol(T, idxs=lastSeg.posM[1]).u
yy = sol(T, idxs=lastSeg.posM[2]).u
zz = sol(T, idxs=lastSeg.posM[3]).u
plt = plot(xx,zz, aspect_ratio=:equal)
display(plt)




end
