using LinearAlgebra
using DifferentialEquations
using SparseArrays
using Plots
using LaTeXStrings
using BenchmarkTools 

N = 3
m = fill(3., N)
#..set spring constant of spring 
k = fill(1., N+1)
for i in 1:N
    k[i] = i
end
l = fill(1., N+1)
D = N+1
#..set damping constant 
c = fill(.2, N)

L = k.*l

#print(L)

#stifness matrix
K = zeros(N, N)

K[1,1] = -k[1] - k[2]
K[1,2] = k[2]
K[N,N] = -k[N] - k[N+1]
K[N,N-1] = k[N]
for i in 2:N-1
    K[i,i] = -k[i] -k[i+1]
    K[i, i-1] = +k[i]
    K[i, i+1] = +k[i+1]
end

print(K)

#..set imposed acceleration on the door in x-direction 
function f(t)
    #return t>=1 
    # return exp(-(t-10)^2/10) 
    return 0.2*sin(t)
end 

#..set imposed acceleration on the door in y-direction 
function g(t)
    #return t>=1 
    return 0 
end
# particular = fill(0., N)

# for i in 1:N
#     particular[i] = (-L[i] + L[i+1])
# end
# particular[N] =  D * k[N+1]
# print(L)
# print(particular)

#..define the right-hand side of the ordinary differential equation of the equation of motion 
function mass_system3!(ddu,du,u,p,t)
    for i in 2:N-1
        ddu[i] = (K[i,i]*u[i] + K[i,i-1]*u[i-1] + K[i,i+1]*u[i+1] - c[i]*du[i] + L[i] - L[i+1])/m[i]
    end

    ddu[1] = (K[1,1]*u[1] + K[1,2]*u[2] - c[1]*du[1] + L[1] - L[2])/m[1] + f(t)
    ddu[N] = (K[N,N]*u[N] + K[N,N-1]*u[N-1] - c[N]*du[N] + L[N] - L[N+1] + D * k[N+1])/m[N]

    # ddu = (K*u -c.*du + particular)./m
    # ddu[1] += f(t)
end

#..set initial position and velocity
u0 = fill(1.,N)
for i in 1:N
    u0[i] = i
end

print(u0)

v0 = zeros(N)
v0[1] = .1
                                    
#..set time begin and end forward
tspan = (0.0,100.0)           

#..define ODE problem to be solved  
prob = SecondOrderODEProblem(mass_system3!,v0,u0,tspan)

#..solve ODE problem 
sol = solve(prob)
#print(sol(1))
#..plot the source term
tvec = Vector(0.:0.01:10.)
fvec = f.(tvec)
# p1 = plot(tvec,fvec,label="Excitation")


#..velocity and position have vars=(1,2) and vars=(3,4), respectively. 
plot(sol,vars=1)
for i in 2:N-1
    plot!(sol,vars=i)
end
p1 = plot!(sol,vars=N)

plot(sol,vars=N+1)
for i in N+2:2*N-1
    plot!(sol,vars=i)
end
p2 = plot!(sol,vars=2*N)


#..plot solution of velocity and position as function of time  
plot(p1,p2,layout=(2,1))
savefig("tuesday0.png")