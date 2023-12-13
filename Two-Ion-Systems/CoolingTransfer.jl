using OrdinaryDiffEq, Plots
using LinearAlgebra

epsilon0 = 8.85*1e-12;
a1 = -0.001
q1 = 0.21801801801801804
Ω_rf = 2*pi*5.2*1e6;
T = 1*1e-3; #1mK
hbar = 1*1e-34;#hbar
kb = 1.380*1e-23; #Boltzmann constant
ω_z1 = omega_rf*sqrt(-a1/2);
ω_r1 = omega_rf/2*sqrt(q1^2/2+a1);
beta = 1/(kb*T);
m1 = 137  * 1.66*1e-27;
m2 = 9000/12 * 1.66*1e-27;
Q1 = 1.6*1e-19;
Q2 = 2*Q1;

mu = m2/m1;
rho = Q2/Q1;

a2 = a1*rho/mu;
q2 = q1*rho/mu;
ω_z2 = Ω_rf*sqrt(-a2/2)
ω_r2 = Ω_rf/2*sqrt(q2^2/2+a2)
println(ω_z1*1e-3/(2*pi),", ",ω_r1*1e-3/(2*pi))
println(ω_z2*1e-3/(2*pi),", ",ω_r2*1e-3/(2*pi))


#!!!!!"1Dimensional problem seems to be working"
# function OneD_Problem(du,u,p,t)
#     z = u[1];
#     dz = u[2];
#     du[1] = dz;
#     du[2] = -ω_z1^2*z
# end
# u0 = [1*1e-6,0]
# E0 = 1/2*m1*ω_z1^2*u0[1]^2
# println("E_0 is ", E0)
# tspan = (0,2*pi/ω_z1)
# prob = ODEProblem(OneD_Problem, u0, tspan)
# sol = solve(prob,Tsit5())
# plot(sol,linewidth = 2,vars = [1,2])
#!!!!!

# #!!!!!! Testing the 2D Problem. Also seems to be working like a charm
# function TwoD_problem(du,u,p,t)
#     z = u[1]
#     dz = u[2]
#     y = u[3]
#     dy = u[4]

#     du[1] = dz
#     du[2] = -ω_z1^2*z
#     du[3] = dy
#     du[4] = -ω_r1^2*y
# end
# tspan = (0,1.5*2*pi/ω_z1);
# prob = ODEProblem(TwoD_problem, u0, tspan)
# u0 = [1*1e-6,0,-1*1e-6,0]
# sol = solve(prob,Tsit5())

# p1 = plot(sol,linewidth = 2 ,vars = [1,3])
# p2 = plot(sol,linewidth = 2, vars = [2,4])
# plot(p1,p2,layout = (2,1),size = (500,500))
# #!!!!!!

# #!!!!! Trying to do simulation of both particles. WITHOUT Coulomb interaction
# function TwoP_NO_COULOMB(du,u,p,t)
#     z1 = u[1]
#     dz1 = u[2]
#     y1 = u[3]
#     dy1 = u[4]
#     z2 = u[5]
#     dz2 = u[6]
#     y2 = u[7]
#     dy2 = u[8]

#     du[1] = dz1
#     du[2] = -ω_z1^2*z1
#     du[3] = dy1
#     du[4] = -ω_r1^2*y1
#     du[5] = dz2
#     du[6] = -ω_z2^2*z2
#     du[7] = dy2
#     du[8] = -ω_r2^2*y2
# end
# v1 = sqrt(2/(m1*beta))
# v2 = -sqrt(2/(m2*beta))
# u0 = [0,v1,0,v1,0,v2,0,v2]
# tspan = (0,1.5*2*pi/ω_z1);
# prob = ODEProblem(TwoP_NO_COULOMB, u0, tspan)
# sol = solve(prob,Tsit5())

# p1 = plot(sol,linewidth = 2,vars = [1,3])
# p2 = plot(sol,linewidth = 2,vars = [2,4])
# p3 = plot(sol,linewidth = 2,vars = [5,7])
# p4 = plot(sol,linewidth = 2,vars = [6,8])
# plot(p1,p2,p3,p4,layout = (2,2))
# #!!!!!

#!!!!!! Adding Coulomb interaction now! (assume z2 < 0, z1 > 0)
z1_eq = cbrt(Q1*Q2/(4*pi*epsilon0*m1*ω_z1^2)/(1+1/rho)^2)
println("z1_eq = ",z1_eq/1e-6," µm")
z2_eq = -z1_eq/rho
println("z2_eq = ",z2_eq/1e-6," µm")
function TwoP_COULOMB(du,u,p,t)
    z1 = u[1]
    dz1 = u[2]
    y1 = u[3]
    dy1 = u[4]
    z2 = u[5]
    dz2 = u[6]
    y2 = u[7]
    dy2 = u[8]
    r = sqrt((z1-z2)^2+(y1-y2)^2)
    r_vec = [z1-z2,y1-y2]/r
    du[1] = dz1
    du[2] = -ω_z1^2*z1+r_vec[1]*Q1*Q2/(4*pi*epsilon0*r^2*m1)
    du[3] = dy1
    du[4] = -ω_r1^2*y1+r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)
    du[5] = dz2
    du[6] = -ω_z2^2*z2-r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
    du[7] = dy2
    du[8] = -ω_r2^2*y2-r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
end
v1 = sqrt(2/(m1*beta))
v2 = -sqrt(2/(m2*beta))
u0 = [z1_eq,v1,0,v1,z2_eq,v2,0,v2]
tspan = (0,3*2*pi/ω_z1);
prob = ODEProblem(TwoP_COULOMB, u0, tspan)
sol = solve(prob,Tsit5(),adaptive = false, dt = 1e-8);

# p1 = plot(sol,linewidth = 2,vars = [1,3])
# p2 = plot(sol,linewidth = 2,vars = [2,4])
# p3 = plot(sol,linewidth = 2,vars = [5,7])
# p4 = plot(sol,linewidth = 2,vars = [6,8])
# plot(p1,p2,p3,p4,layout = (2,2))

p1 = plot(sol.t,[sol[1,:].-z1_eq,sol[5,:].-z2_eq],linewidth = 2)
p2 = plot(sol.t,[sol[3,:],sol[7,:]],linewidth = 2)
plot(p1,p2,layout = (2,1))
