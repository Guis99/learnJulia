include("fe.jl")

using Plots

function main()
    numElem = 4
    order = 4

    gs = 0.0
    ge = 2.0

    allSize = (ge - gs) / numElem

    elemSizes = fill(allSize,numElem)
    mesh = makeMesh(numElem,4,gs,ge, elemSizes)

    nNodes = length(mesh.Nodes)[1]

    LHSCoeff = 105.0
    RHSCoeff = 1.0

    bc = boundaryCondition1D([false,true], [0,10])
    display(bc)

    sol = solveProblem(numElem, order, gs, ge, LHSCoeff, RHSCoeff, elemSizes, bc)

    display(sol)

    plot(getNodePos(mesh, collect(1:nNodes)), sol)
end

main()
