using LinearAlgebra
using DifferentialEquations
using SparseArrays
using Plots
using LaTeXStrings
using BenchmarkTools 

#..set mass of point mass 
m = 5 
#..set spring constant of spring 
k = 1.0 
#..set damping constant 
c = 0.5 

#..set imposed acceleration on the door 
function f(t)
    #return exp(-(t-2)^2/0.01)
    return t>=2 
end 

#..define the right-hand side of the ordinary differential equation of the equation of motion 
function mass_system(du,u,p,t)
    # solve \ddot{u} = -(k/m) u - (c/m) \dot u + f(t) 
    ddu = -(k/m)*u - (c/m)*du + f(t)
end

function mass_system!(ddu,du,u,p,t)
    # solve \ddot{u} = -(k/m) u - (c/m) \dot u + f(t) 
    ddu = -(k/m)*u - (c/m)*du + f(t)
end

#..set initial position and velocity
u0 = 1.0                                      
v0 = 0.0 
#..set time begin and end forward
tspan = (0.0,10.0)               

#..define ODE problem to be solved  
prob = SecondOrderODEProblem(mass_system,v0,u0,tspan)

#..solve ODE problem 
sol = solve(prob)

#..plot the source term
tvec = Vector(0.:0.01:10.)
fvec = f.(tvec)
p1 = plot(tvec,fvec,label="Excitation")

#..plot solution of velocity and position as function of time  
plot(sol,vars=1,label="Velocity")
p2 = plot!(sol,vars=2,title="m = 5, k = 1.0, c = 0.5",label="Position")

plot(p1,p2,layout=(2,1))
print("gonna save")
savefig("monday1.png")
print("end")