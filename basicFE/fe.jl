using Polynomials
using LinearAlgebra

struct Node
    position::Float64
end

struct Element1D
    leftBoundary::Float64
    rightBoundary::Float64
    nodeIdxs::Vector{Int}
end

struct Mesh1D
    numElem::Int
    order::Int
    domainStart::Float64
    domainEnd::Float64
    polynomialFuncs::Vector{Polynomial}
    Elements::Vector{Element1D}
    elemSizes::Vector{Float64}
    Nodes::Vector{Node}
end

struct boundaryCondition1D
    isDir::Vector{Bool}
    dirVal::Vector{Float64}
end

function getNodePos(mesh::Mesh1D, nodeIdxs::Vector{Int64})
    n = length(nodeIdxs)[1]
    out = Vector{Float64}(undef, n)
    for i=1:n
        out[i] = mesh.Nodes[nodeIdxs[i]].position
    end

    return out
end

function genGaussPoints(order)
    numPoints = order + 1
    evalPts = [pi * range(0,1,numPoints);]
    gaussPoints = reverse(cos.(evalPts))
    return gaussPoints
end

function genLagrange(gaussPoints::Vector{Float64})
    numPts = size(gaussPoints)[1]
    out = Vector{Polynomial}(undef, numPts)
    for i in 1:numPts
        excluded = [gaussPoints[1:i-1]; gaussPoints[i+1:end]]
        poly = fromroots(excluded)
        out[i] = poly / poly(gaussPoints[i])
    end
    return out
end

function firstDerivatives(basisFuncs::Vector{Polynomial})
    return derivative.(basisFuncs)
end

function firstDerivativesMat(dbasisFuncsdx::Vector{Polynomial{Float64, :x}}, gaussPoints::Vector{Float64})
    n = length(dbasisFuncsdx)[1]
    out = Array{Float64}(undef, (n,n))

    # for i=1:n
    #     deriv = dbasisFuncsdx[i]
    #     evalD = deriv.(gaussPoints)
    #     out[:,i] = evalD
    # end

    foreach(i -> out[:,i] = dbasisFuncsdx[i].(gaussPoints), 1:n |> Array)

    return out
end

function getWeights(basisFuncs::Vector{Polynomial})
    numBasis = size(basisFuncs)[1]
    out = Vector{Float64}(undef,numBasis)
    ADs = integrate.(basisFuncs)

    for i=1:numBasis
        currBasis = ADs[i]
        out[i] = currBasis(1) - currBasis(-1)
    end

    return out
end

function doTransform(gaussPoints, a::Float64, b::Float64)
    ao = -1
    bo = 1

    intlength = b - a
    intlengtho = bo - ao
    newCoords = @. a + (gaussPoints - ao) * intlength / intlengtho

    return newCoords
end

function makeMesh(numElem, order, start::Float64, endpoint::Float64, elemSizes::Vector{Float64})
    gaussPoints = genGaussPoints(order)
    basisfuncs = genLagrange(gaussPoints)
    nNodes = numElem * order + 1

    nodes = Vector{Node}(undef, nNodes)
    elems = Vector{Element1D}(undef, numElem)

    nodes[1] = Node(start)

    currPos = start

    for i=0:numElem-1
        elemSize = elemSizes[i+1]
        idxStart = order * i + 1
        nodeIdxs = idxStart:idxStart+order |> Array
        elems[i+1] = Element1D(currPos, currPos+elemSize, nodeIdxs)

        gaussPointsTransformed = doTransform(gaussPoints, currPos, currPos+elemSize)
        for j=2:order+1
            globalIdx = idxStart + j - 1
            nodes[globalIdx] = Node(gaussPointsTransformed[j])
        end
        currPos += elemSize
    end

    mesh = Mesh1D(numElem, order, start, endpoint, basisfuncs, elems, elemSizes, nodes)
    return mesh
end

function assembleStiffnessMatrix(mesh::Mesh1D, LHSCoeff::Function, RHSCoeff::Function)
    n = length(mesh.Nodes)[1]
    outK = zeros(n,n)
    outF = zeros(n,1)

    display(outK)

    gaussPoints = genGaussPoints(mesh.order)

    firstDerivs = firstDerivatives(mesh.polynomialFuncs)
    firstDerivsMat = firstDerivativesMat(firstDerivs, gaussPoints)
    weights = getWeights(mesh.polynomialFuncs)
    weightsMat = diagm(weights)

    for elem::Element1D=mesh.Elements
        jac = 2 / (elem.rightBoundary - elem.leftBoundary) # dxdz

        nodeIdxs = elem.nodeIdxs
        nodePos = getNodePos(mesh, nodeIdxs)

        localElemK = jac .* firstDerivsMat' * weightsMat * diagm(LHSCoeff(nodePos)) * firstDerivsMat
        localElemF = weightsMat * RHSCoeff(nodePos) ./ jac

        outK[nodeIdxs,nodeIdxs] += localElemK
        outF[nodeIdxs] += localElemF
    end

    return (outK, outF)
end

function solveProblem(numElem, order, domainStart, domainEnd, LHSCoeff, RHSCoeff, elemSizes, boundaryConds::boundaryCondition1D)
    n = numElem * order + 1
    if !isa(LHSCoeff, Function)
        LHSCoeff_f = x -> fill(LHSCoeff, length(x)[1])
    else
        LHSCoeff_f = LHSCoeff
    end

    if !isa(RHSCoeff, Function)
        RHSCoeff_f = x -> fill(RHSCoeff, length(x)[1])
    else
        RHSCoeff_f = RHSCoeff
    end

    mesh = makeMesh(numElem, order, domainStart, domainEnd, elemSizes)
    (K, f) = assembleStiffnessMatrix(mesh, LHSCoeff_f, RHSCoeff_f)

    isDirVec = boundaryConds.isDir
    numDirBC = sum(isDirVec)

    idxsToRemove = [1,n]

    idxsToRemove = idxsToRemove[isDirVec]
    idxsToKeep = setdiff(collect(1:n),idxsToRemove)

    display(idxsToKeep)
    display(idxsToRemove)

    CS = zeros(n,n-numDirBC)
    NS = zeros(n,numDirBC)
    
    for i=1:length(idxsToRemove)[1]
        NS[idxsToRemove[i],i] = 1.0
    end

    for i=1:length(idxsToKeep)[1]
        CS[idxsToKeep[i],i] = 1.0
    end

    K11 = CS' * K * CS
    K12 = CS' * K * NS
    F11 = CS' * f

    boundaryVals = boundaryConds.dirVal[isDirVec]

    solInterm =  K11 \ (F11 - K12 * boundaryVals)
    sol = CS * solInterm + NS * boundaryVals

    return sol
end