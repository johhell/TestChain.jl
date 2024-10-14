"""
einzelne Segmente, die mauall verbunden werden
    function addSegments!(kette, settings)
"""
module Kette

using Modia
using StaticArrays
using LinearAlgebra
using DiffRules: @define_diffrule
using Plots



include("Tparameters.jl")

vec3Null = SVector(0.0,0.0,0.0)


se = Settings(;
        segments=5,
        v_ro = 0.0,
        frictionCoeff = 0.0,
        endForce=vec3Null,
    )


Pin3 = Model(
    position = Var(potential=true, start= vec3Null),
    force = Var(flow=true, start= vec3Null),
)


Fixed_P0   = Model( a0 = Pin3, equations = :[ a0.position = [0.0, 0.0,0.0]] ) # kein SVector !! sonst Problem bei diff

Ende_F0    = Model(
    endForce = Par(value = se.endForce),
    a0 = Pin3,
    equations = :[ a0.force = endForce ],
)


Distance(a) = norm(a)
DDistance(a) = [a[1] a[2] a[3]] / norm(a)
@define_diffrule Base.Distance(a)  = :( DDistance($a) )


Segment = Model(
    mass        = se.mass_per_seg,
    length0     = se.len_per_segment,
    cSpring     = se.c_spring0,
    dSpring     = se.damping,
    a1 = Pin3,
    a2 = Pin3,
    gravity     = se.g_earth,
    posM        = Var(init= vec3Null),
    velM        = Var(init= vec3Null),
    accM        = Var(start=vec3Null),
    forceSpring = Var(start=vec3Null),
    segment     = Var(start=vec3Null),
    forceSpringABS  = Var(),
    velABS = Var(),
    segABS = Var(),

    equations = :[
        segment = a1.position-posM
        forceSpringABS = (segABS-length0)*cSpring + velABS*dSpring
        segABS = Distance(segment)
        velABS = (der(segABS))[1]   # Ergebnis ist Vector mit 1 Element
        segUnit = segment/segABS        # Richtung Kraft - Länge=1.0
        velM = der(posM)
        accM = der(velM)
        forceSpring = segUnit*forceSpringABS
        (accM-gravity)*mass = forceSpring + a2.force
        a1.force = forceSpring
        a2.position = posM
    ],
)


kette = Model(
    templateSeg = Segment, # dummy
    fix = Fixed_P0,
    ende = Ende_F0,
    equations = :[   ],
    connect = :[ ]
)


function addSegments!(kette, settings)
    liste = Vector{String}()
    for i = 1:settings.segments
        name = "Seg$i"
        segSym = Symbol(name)
        push!(liste,name)
        # neues Segment
        kette[segSym] = deepcopy(kette[:templateSeg])
        # Startwerte POSITION
        kette[segSym][:posM][:init] = settings.positionStart[i]
    end
    # Platzhalter entfernt
    delete!(kette,:templateSeg)

    function addConn(txt)
        push!(conns, Meta.parse(txt))
    end

    # CONNECT
    # neue Liste
    conns = kette[:connect].args=Vector{Any}()
    # Anfang & Ende
    addConn("fix.a0, $(liste[1]).a1")
    addConn("ende.a0, $(liste[end]).a2")
    # Verbindungen Segmente
    for i = 2:settings.segments
        addConn("$(liste[i-1]).a2, $(liste[i]).a1")
    end
    liste
end


segListe = addSegments!(kette, se)


instModel = @instantiateModel(kette,
    unitless=true,
    log=false,
    logTiming=false,
    logDetails=false,
    logCode=false,
    logStateSelection=false,
#     saveCodeOnFile="CODE.txt.jl",
    )


function Simulation(inst)
    simulate!(inst,
        stopTime = se.duration,
        interval = 0.02,
        log = false,
        logStates=true,
        tolerance = 1e-6  )
end


@time res = Simulation(instModel)


function Bilder(segs, zeit)
    p1 = plot(;xlabel="x", ylabel="z")
    p2 = plot(;xlabel="time", ylabel="x")
    p3 = plot(;xlabel="time", ylabel="z")
    p4 = plot(;xlabel="time", ylabel="Δ length")
    for b in segs
        posM = get_result(instModel,"$(b).posM")
        segABS = get_result(instModel,"$(b).segABS") .- se.len_per_segment
        force = get_result(instModel,"$(b).forceSpringABS")
        plot!(p1,posM[:,1],posM[:,3], label="")
        scatter!(p1,[posM[end,1]],[posM[end,3]], label="")
        plot!(p2,zeit, posM[:,1], label="")
        plot!(p3,zeit, posM[:,3], label="")
        plot!(p4,zeit, segABS, label="")
    end
    plt=plot(p1,p2,p3,p4)
end

# plt = Bilder(segListe, res.t)
# display(plt)

# Theorie für 1 Segment
# @show T0 = 2pi*sqrt(se.l0/9.81)

if false
    FH = open("EndPosition.txt","a")
    b = segListe[end]
    posM = get_result(instModel,"$(b).posM")
    write(FH,"$(se.segments)  $(posM[end,1])   $(posM[end,2])   $(posM[end,3])  \n"  )
    close(FH)
end



# position of the last particle
T = res.t


using HDF5

h5open("SegModiaL0.hdf5", "w") do fid
    fid["zeit"] = T
    fid["model"] = "Seg-Modia0"
    fid["segments"] = se.segments
    gruppe= create_group(fid, "positions")
    for i = 1:se.segments
        b = segListe[i]
        posM = get_result(instModel,"$(b).posM")
#         xx = posM[:,1]
#         yy = posM[:,2]
#         zz = posM[:,3]
#         gruppe["pos$(i)"] = hcat(xx,yy,zz)
        gruppe["pos$(i)"] = posM    # xx,yy,zz in Spalten

    end
end


using Plots

b = segListe[end]
posM = get_result(instModel,"$(b).posM")
xx = posM[:,1]
yy = posM[:,2]
zz = posM[:,3]

plt = plot(xx,zz, aspect_ratio=:equal)
display(plt)



end
