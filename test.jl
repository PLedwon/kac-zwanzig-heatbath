using DifferentialEquations
println("done")
alpha = 0.5 #Setting alpha to 1/2
f(y,p,t) = alpha*y
u0 = 1.5
timespan = (0.0,1.0) # Solve from time = 0 to time = 1
println("done")
prob = ODEProblem(f,u0,timespan)
println("done")
#sol = solve(prob) # Solves the ODE
#plot(sol) # Plots the solution using Plots.jl
