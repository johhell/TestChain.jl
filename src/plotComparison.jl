
module bilder

using Plots
using Glob


using TestChain


liste = glob("*.hdf5")

dataList = Vector{Daten}()
for f in liste
    d = Daten(f)
    push!(dataList,d)
end



function BildXZ(Dlist::Vector{Daten}, plt)    # nur 1 Element
    for b in Dlist
        P1 = b.positions[end]
        plot!(plt, P1.xx, P1.zz, label=b.name)
    end
end


function BildX(Dlist::Vector{Daten}, plt)    # nur 1 Element
    for b in Dlist
        P1 = b.positions[end]
        plot!(plt, b.time, P1.xx, label=b.name)
    end
end

function BildZ(Dlist::Vector{Daten}, plt)    # nur 1 Element
    for b in Dlist
        P1 = b.positions[end]
        plot!(plt, b.time, P1.zz, label=b.name)
    end
end



p1 = plot(;xlabel="x", ylabel="z")
p2 = plot(;xlabel="time", ylabel="x")
p3 = plot(;xlabel="time", ylabel="z")



BildXZ(dataList, p1)
BildX(dataList, p2)
BildZ(dataList, p3)


p = plot(p1,p2,p3, size=(800,500), plot_title = "comparision level 0")


pngName = "../images/CompLevel0.png"
printstyled("savePNG: $(pngName)\n", color= :blue)
savefig(p, pngName)
display(p)

end
