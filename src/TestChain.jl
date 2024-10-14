module TestChain

using HDF5

export ElementPos, Daten


struct ElementPos
    nummer::Int64
    xx::Vector{Float64}
    yy::Vector{Float64}
    zz::Vector{Float64}
    function ElementPos(nr::Int64, xyz::Matrix{Float64})
        x = xyz[:,1]
        y = xyz[:,2]
        z = xyz[:,3]
        new(nr, x,y,z)
    end
end


struct Daten
    name::String
    time ::Vector{Float64}
    segments::Int64
    positions::Vector{ElementPos}
    function Daten(hdf)
        FH = h5open("$(hdf)", "r")
        t = read(FH, "zeit")
        n = read(FH, "model")
        s = read(FH, "segments")
        p = Vector{ElementPos}()
        gruppe= read(FH, "positions")
        for i = 1:s
            pxyz = gruppe["pos$(i)"]
            push!(p,ElementPos(i,pxyz))
        end
        new(n, t, s, p)
    end
end


end
