using Polynomials
include("fe.jl")

p=fromroots([1, 2, 3])
gp1=genGaussPoints(1)
gp5=genGaussPoints(5)
gp10=genGaussPoints(10)
display(p)

pl1 = fromroots(gp1)
pl5 = fromroots(gp5)
pl10 = fromroots(gp10)

display(pl1)
display(pl1 / 10)
display(pl5)
display(pl10)

display(pl1.([-1,0,1]))

display([pl1 pl5 pl10])


dpl1dx=derivative(pl1)
display(dpl1dx)