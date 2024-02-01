using OrdinaryDiffEq, Plots
using LinearAlgebra
using LaTeXStrings
using DifferentialEquations
using Optim


epsilon0 = 8.85*1e-12;
a1 = -0.001
q1 = 0.21801801801801804
omega_rf = 2*pi*5.2*1e6;
T = 2*1e-3; #1mK
hbar = 1*1e-34;#hbar
kb = 1.380*1e-23; #Boltzmann constant
ω_z1 = omega_rf*sqrt(-a1/2);
ω_r1 = omega_rf/2*sqrt(q1^2/2+a1);
beta = 1/(kb*T);
m1 = 137  * 1.66*1e-27;
m2 = 9000 * 1.66*1e-27;
Q1 = 1.6*1e-19;
Q2 = 24*Q1;
V_DC = -a1*m1*(3/4*2.7*1e-3)^2*omega_rf^2/(4*Q1*0.248)
println("DC Voltage is:    ", V_DC)
z0 = 3/4*2.7*1e-3
r0 = 3/4*3.5*1e-3

mu = m2/m1;
rho = Q2/Q1;

a2 = a1*rho/mu;
q2 = q1*rho/mu;
ω_z2 = omega_rf*sqrt(-a2/2)
ω_r2 = omega_rf/2*sqrt(q2^2/2+a2)
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

# #!!!!!! Adding Coulomb interaction now! (assume z2 < 0, z1 > 0)
# z1_eq = cbrt(Q1*Q2/(4*pi*epsilon0*m1*ω_z1^2)/(1+1/rho)^2)
# println("z1_eq = ",z1_eq/1e-6," µm")
# z2_eq = -z1_eq/rho
# println("z2_eq = ",z2_eq/1e-6," µm")
# function TwoP_COULOMB(du,u,p,t)
#     z1 = u[1]
#     dz1 = u[2]
#     y1 = u[3]
#     dy1 = u[4]
#     z2 = u[5]
#     dz2 = u[6]
#     y2 = u[7]
#     dy2 = u[8]
#     r = sqrt((z1-z2)^2+(y1-y2)^2)
#     r_vec = [z1-z2,y1-y2]/r
#     du[1] = dz1
#     du[2] = -ω_z1^2*z1+r_vec[1]*Q1*Q2/(4*pi*epsilon0*r^2*m1)
#     du[3] = dy1
#     du[4] = -ω_r1^2*y1+r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)
#     du[5] = dz2
#     du[6] = -ω_z2^2*z2-r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
#     du[7] = dy2
#     du[8] = -ω_r2^2*y2-r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
# end
# v1 = sqrt(2/(m1*beta))
# v2 = -sqrt(2/(m2*beta))
# u0 = [z1_eq,v1,0,v1,z2_eq,v2,0,v2]
# tspan = (0,3*2*pi/ω_z1);
# prob = ODEProblem(TwoP_COULOMB, u0, tspan)
# sol = solve(prob,Tsit5(),adaptive = false, dt = 1e-8);

# # p1 = plot(sol,linewidth = 2,vars = [1,3])
# # p2 = plot(sol,linewidth = 2,vars = [2,4])
# # p3 = plot(sol,linewidth = 2,vars = [5,7])
# # p4 = plot(sol,linewidth = 2,vars = [6,8])
# # plot(p1,p2,p3,p4,layout = (2,2))


# K11 = ω_z1^2*(1+2/(1+1/rho)) 
# K22 = ω_z2^2*(1+2/(1+rho))
# K12 = -2*ω_z1^2/(sqrt(mu)*(1+1/rho))
# z_MAT = [K11 K12 ; K12 K22]
# egz = eigen(z_MAT);

# K33 = ω_r1^2-ω_z1^2/(1+1/rho);
# K44 = ω_r2^2-ω_z1^2/(mu*(1+1/rho));
# K34 = -1/2*K12

# r_MAT = [K33 K34; K34 K44];
# egy = eigen(r_MAT);
# if sign(egy.vectors[1,1]) == sign(egy.vectors[1,2])
#     ω_yin = sqrt(egy.values[1]);
#     ω_yout = sqrt(egy.values[2]);
#     vector_yin = egy.vectors[1,:];
#     vector_yout = egy.vectors[2,:];
# else
#     ω_yin = sqrt(egy.values[2]);
#     ω_yout = sqrt(egy.values[1]);
#     vector_yin = egy.vectors[2,:];
#     vector_yout = egy.vectors[1,:];
# end
# if sign(egz.vectors[1,1]) == sign(egz.vectors[1,2])
#     ω_zin = sqrt(egz.values[1]);
#     ω_zout = sqrt(egz.values[2]);
#     vector_z_in = egz.vectors[1,:];
#     vector_z_out = egz.vectors[2,:];
# else
#     ω_zin = sqrt(egz.values[2]);
#     ω_zout = sqrt(egz.values[1]);
#     vector_z_in = egz.vectors[2,:];
#     vector_z_out = egz.vectors[1,:];
# end

# chi_in = zeros(1,length(sol[1,:]));
# chi_out = zeros(1,length(sol[1,:]));
# ts = zeros(1,length(sol.t))
# for j in eachindex(chi_in)
#     ts[j] = sol.t[j]
#     z1 = sol[1,j]
#     z2 = sol[5,j]
#     chi_in[j] = (z1-z1_eq)*vector_z_in[1]+(z2-z2_eq)*vector_z_in[2];
#     chi_out[j] = (z1-z1_eq)*vector_z_out[1]+(z2-z2_eq)*vector_z_out[2];
# end
# p1 = plot(sol.t,[sqrt(137).*(sol[1,:].-z1_eq).*vector_z_in[1]+sqrt(9000/12).*(sol[5,:].-z2_eq).*vector_z_in[2],sqrt(137).*(sol[1,:].-z1_eq)*vector_z_out[1]+sqrt(9000/12).*(sol[5,:].-z2_eq).*vector_z_out[2]],linewidth = 2)
# p2 = plot(sol.t,[sqrt(137).*sol[3,:].*vector_yin[1]+sqrt(9000/12).*sol[7,:].*vector_yin[2],sqrt(137).*sol[3,:].*vector_yout[1]+sqrt(9000/12).*sol[7,:].*vector_yout[2]],linewidth = 2)
# plot(p1,p2,layout = (2,1))
# #!!!!!!!!!!!!!




# #!!!!!!!!!!! code with doppler cooling as well
# z1_eq = cbrt(Q1*Q2/(4*pi*epsilon0*m1*ω_z1^2)/(1+1/rho)^2)
# println("z1_eq = ",z1_eq/1e-6," µm")
# z2_eq = -z1_eq/rho
# println("z2_eq = ",z2_eq/1e-6," µm")
# Sat_ratio = 0.1
# k = 12.7 *1e6
# gamma = 2*pi*15.2*1e6
# delta = -gamma/2
# cooling_friction = 4*hbar*k^2 * Sat_ratio*(-2*delta/gamma)/(1+(2*delta/gamma)^2)^2
# println("Cooling friction factor is:   ",cooling_friction)

# function TwoP_COULOMB_DOPPL(du,u,p,t)
#     z1 = u[1]
#     dz1 = u[2]
#     y1 = u[3]
#     dy1 = u[4]
#     z2 = u[5]
#     dz2 = u[6]
#     y2 = u[7]
#     dy2 = u[8]
#     r = sqrt((z1-z2)^2+(y1-y2)^2)
#     r_vec = [z1-z2,y1-y2]/r
#     du[1] = dz1
#     du[2] = -ω_z1^2*z1+r_vec[1]*Q1*Q2/(4*pi*epsilon0*r^2*m1)-1/sqrt(2)*cooling_friction*u[2]/m1
#     du[3] = dy1
#     du[4] = -ω_r1^2*y1+r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)-1/sqrt(2)*cooling_friction*u[4]/m1
#     du[5] = dz2
#     du[6] = -ω_z2^2*z2-r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
#     du[7] = dy2
#     du[8] = -ω_r2^2*y2-r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
# end
# v1 = sqrt(2/(m1*beta))
# v2 = -sqrt(2/(m2*beta))
# u0 = [z1_eq,v1,0,v1,z2_eq,v2,0,v2]
# tspan = (0,50*2*pi/ω_z1);
# prob = ODEProblem(TwoP_COULOMB_DOPPL, u0, tspan)
# sol = solve(prob,Tsit5(),adaptive = false, dt = 1e-8);

# # p1 = plot(sol,linewidth = 2,vars = [1,3])
# # p2 = plot(sol,linewidth = 2,vars = [2,4])
# # p3 = plot(sol,linewidth = 2,vars = [5,7])
# # p4 = plot(sol,linewidth = 2,vars = [6,8])
# # plot(p1,p2,p3,p4,layout = (2,2))


# K11 = ω_z1^2*(1+2/(1+1/rho)) 
# K22 = ω_z2^2*(1+2/(1+rho))
# K12 = -2*ω_z1^2/(sqrt(mu)*(1+1/rho))
# z_MAT = [K11 K12 ; K12 K22]
# egz = eigen(z_MAT);

# K33 = ω_r1^2-ω_z1^2/(1+1/rho);
# K44 = ω_r2^2-ω_z1^2/(mu*(1+1/rho));
# K34 = -1/2*K12

# r_MAT = [K33 K34; K34 K44];
# egy = eigen(r_MAT);
# if sign(egy.vectors[1,1]) == sign(egy.vectors[1,2])
#     ω_yin = sqrt(egy.values[1]);
#     ω_yout = sqrt(egy.values[2]);
#     vector_yin = egy.vectors[1,:];
#     vector_yout = egy.vectors[2,:];
# else
#     ω_yin = sqrt(egy.values[2]);
#     ω_yout = sqrt(egy.values[1]);
#     vector_yin = egy.vectors[2,:];
#     vector_yout = egy.vectors[1,:];
# end
# if sign(egz.vectors[1,1]) == sign(egz.vectors[1,2])
#     ω_zin = sqrt(egz.values[1]);
#     ω_zout = sqrt(egz.values[2]);
#     vector_z_in = egz.vectors[1,:];
#     vector_z_out = egz.vectors[2,:];
# else
#     ω_zin = sqrt(egz.values[2]);
#     ω_zout = sqrt(egz.values[1]);
#     vector_z_in = egz.vectors[2,:];
#     vector_z_out = egz.vectors[1,:];
# end
# # chi_in = zeros(1,length(sol[1,:]));
# # chi_out = zeros(1,length(sol[1,:]));
# # ts = zeros(1,length(sol.t))
# # for j in eachindex(chi_in)
# #     ts[j] = sol.t[j]
# #     z1 = sol[1,j]
# #     z2 = sol[5,j]
# #     chi_in[j] = (z1-z1_eq)*vector_z_in[1]+(z2-z2_eq)*vector_z_in[2];
# #     chi_out[j] = (z1-z1_eq)*vector_z_out[1]+(z2-z2_eq)*vector_z_out[2];
# # end
# # p1 = plot(sol.t,[sqrt(137).*(sol[1,:].-z1_eq).*vector_z_in[1]+sqrt(9000/12).*(sol[5,:].-z2_eq).*vector_z_in[2],sqrt(137).*(sol[1,:].-z1_eq)*vector_z_out[1]+sqrt(9000/12).*(sol[5,:].-z2_eq).*vector_z_out[2]],linewidth = 2)
# # p2 = plot(sol.t,[sqrt(137).*sol[3,:].*vector_yin[1]+sqrt(9000/12).*sol[7,:].*vector_yin[2],sqrt(137).*sol[3,:].*vector_yout[1]+sqrt(9000/12).*sol[7,:].*vector_yout[2]],linewidth = 2)
# # plot(p1,p2,layout = (2,1))

# #Chi is axial, xi is radial
# chi_v_in = sqrt(m1).*(sol[2,:].*vector_z_in[1]).+sqrt(m2).*sol[6,:].*vector_z_in[2];
# chi_v_out = sqrt(m1).*sol[2,:].*vector_z_out[1].+sqrt(m2).*sol[6,:].*vector_z_out[2];
# xi_v_in = sqrt(m1).*sol[4,:].*vector_yin[1].+sqrt(m2).*sol[8,:].*vector_yin[2];
# xi_v_out = sqrt(m1).*sol[4,:].*vector_yout[1].+sqrt(m2).*sol[8,:].*vector_yout[2];



# Tz_in = 1/2 .*chi_v_in.^2;
# Tz_out = 1/2 .*chi_v_out.^2;
# Tr_in = 1/2 .*xi_v_in.^2
# Tr_out = 1/2 .*xi_v_out.^2
# T1 = 1/2 * m1 .* sol[2,:].^2
# T2 = 1/2 * m2 .* sol[6,:].^2
# p1 = plot(sol.t,[Tz_in./(kb*T)], label = L"T_{in,z}")
# p2 = plot(sol.t,[Tz_out./(kb*T)],label = L"T_{out,z}")
# p3 = plot(sol.t,[Tr_in./(kb*T)],label = L"T_{in,r}")
# p4 = plot(sol.t,[Tr_out./(kb*T)],label = L"T_{out,r}")
# plot(p1,p2,p3,p4,layout = (2,2))
# #!!!!!!!!!!!


#!!!!!!!!!!! code with doppler cooling + attempted transfer
z1_eq = cbrt(Q1*Q2/(4*pi*epsilon0*m1*ω_z1^2)/(1+1/rho)^2)
println("z1_eq = ",z1_eq/1e-6," µm")
z2_eq = -z1_eq/rho
println("z2_eq = ",z2_eq/1e-6," µm")
Sat_ratio = 0.1
k = 12.7 *1e6
gamma = 2*pi*15.2*1e6
delta = -gamma/2
cooling_friction = 4*hbar*k^2 * Sat_ratio*(-2*delta/gamma)/(1+(2*delta/gamma)^2)^2
println("Cooling friction factor is:   ",cooling_friction)
eta0 = sqrt(z0^2+r0^2);
V0 = 1/10;
K11 = ω_z1^2*(1+2/(1+1/rho)) 
K22 = ω_z2^2*(1+2/(1+rho))
K12 = -2*ω_z1^2/(sqrt(mu)*(1+1/rho))
z_MAT = [K11 K12 ; K12 K22]
egz = eigen(z_MAT);

K33 = ω_r1^2-ω_z1^2/(1+1/rho);
K44 = ω_r2^2-ω_z1^2/(mu*(1+1/rho));
K34 = -1/2*K12
r_MAT = [K33 K34; K34 K44];
egy = eigen(r_MAT);
if sign(egy.vectors[1,1]) == sign(egy.vectors[1,2])
    ω_yin = sqrt(egy.values[1]);
    ω_yout = sqrt(egy.values[2]);
    vector_yin = egy.vectors[1,:];
    vector_yout = egy.vectors[2,:];
else
    ω_yin = sqrt(egy.values[2]);
    ω_yout = sqrt(egy.values[1]);
    vector_yin = egy.vectors[2,:];
    vector_yout = egy.vectors[1,:];
end
if sign(egz.vectors[1,1]) == sign(egz.vectors[1,2])
    ω_zin = sqrt(egz.values[1]);
    ω_zout = sqrt(egz.values[2]);
    vector_z_in = egz.vectors[1,:];
    vector_z_out = egz.vectors[2,:];
else
    ω_zin = sqrt(egz.values[2]);
    ω_zout = sqrt(egz.values[1]);
    vector_z_in = egz.vectors[2,:];
    vector_z_out = egz.vectors[1,:];
end

function TwoP_COULOMB_DOPPL(du,u,p,t)
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
    du[2] = -ω_z1^2*z1+r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)-1/sqrt(2)*cooling_friction*u[2]/m1#+V0*(z1/eta0^2*Q1)/m1*sin((ω_zin-ω_zout)*t)
    du[3] = dy1
    du[4] = -ω_r1^2*y1+r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)-1/sqrt(2)*cooling_friction*u[4]/m1
    du[5] = dz2
    du[6] = -ω_z2^2*z2-r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)#+V0*(z2/eta0^2*Q2)/m2*sin((ω_zin-ω_zout)*t)
    du[7] = dy2
    du[8] = -ω_r2^2*y2-r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
end
v1 = sqrt(2/(m1*beta))
v2 = -sqrt(2/(m2*beta))
#u0 = [z1_eq+1e-18/sqrt(m1)*vector_z_in[1],0,0,0,z2_eq+1e-18/sqrt(m2)*vector_z_in[2],0,0,0]
u0 = [z1_eq+1e-18/sqrt(m1)*vector_z_out[1],0,0,0,z2_eq+1e-18/sqrt(m2)*vector_z_out[2],0,0,0]

tspan = (0,1000*2*pi/ω_zin);
prob = ODEProblem(TwoP_COULOMB_DOPPL, u0, tspan)
sol = solve(prob,Tsit5(),dtmax = 1e-7);
#Transfer rate at these settings is approx 300Hz, this means period of approx. 20ms.
# p1 = plot(sol,linewidth = 2,vars = [1,3])
# p2 = plot(sol,linewidth = 2,vars = [2,4])
# p3 = plot(sol,linewidth = 2,vars = [5,7])
# p4 = plot(sol,linewidth = 2,vars = [6,8])
# plot(p1,p2,p3,p4,layout = (2,2))
println(vector_z_in)
#Chi is axial, xi is radial
chi_v_in = sqrt(m1).*(sol[2,:].*vector_z_in[1]).+sqrt(m2).*sol[6,:].*vector_z_in[2];
chi_v_out = sqrt(m1).*sol[2,:].*vector_z_out[1].+sqrt(m2).*sol[6,:].*vector_z_out[2];
xi_v_in = sqrt(m1).*sol[4,:].*vector_yin[1].+sqrt(m2).*sol[8,:].*vector_yin[2];
xi_v_out = sqrt(m1).*sol[4,:].*vector_yout[1].+sqrt(m2).*sol[8,:].*vector_yout[2];
U0 = Q1*Q2/(4*pi*epsilon0*(z1_eq-z2_eq))+1/2*m1*ω_z1^2*z1_eq^2+1/2*m2*ω_z2^2*z2_eq^2
println(U0/(kb*T))
# plot1 = plot(sol.t,chi_v_in)
# plot2 = plot(sol.t,chi_v_out)

p1 = plot(sol.t,[sqrt(137).*(sol[1,:].-z1_eq).*vector_z_in[1]+sqrt(9000).*(sol[5,:].-z2_eq).*vector_z_in[2],sqrt(137).*(sol[1,:].-z1_eq)*vector_z_out[1]+sqrt(9000).*(sol[5,:].-z2_eq).*vector_z_out[2]],linewidth = 2)



println(v1*k/(gamma))
T_1 = 1/2 * m1 .*sol[2,:].^2
T_2 = 1/2 * m2 .*sol[6,:].^2
T_3 = 1/2 *m1 .*sol[4,:].^2
T4 = 1/2 *m2 .*sol[8,:].^2
U1 = 1/2* m1 *ω_z1^2 .*sol[1,:].^2
U2 = 1/2* m2 *ω_z2^2 .*sol[5,:].^2
U3 = 1/2* m1 *ω_r1^2 .*sol[3,:].^2
U4 = 1/2* m2 *ω_r2^2 .*sol[7,:].^2
println(ω_zin/2/pi*1e-3)
println(ω_zout/2/pi*1e-3)
UC = Q1*Q2./(4*pi*epsilon0.*sqrt.((sol[1,:].-sol[5,:]).^2+(sol[3,:].-sol[7,:]).^2))
# Ut = sol[1,:].^2 .*V0*(1/eta0^2*Q1)/m1*sin((ω_zin-ω_zout).*sol.t)+sol[1,:].^2 .*V0*(1/eta0^2*Q1)/m1*sin((ω_zin-ω_zout).*sol.t)
 #Checked that energy is conserved, when no damping is present
#!!!!!!!!!



function OptimizerFunc(delta)
    function opt_helper(du,u,p,t)
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
        du[2] = -ω_z1^2*z1+r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)+V0*(z1/eta0^2*Q1)/m1*sin((ω_zin-ω_zout+delta)*t)
        du[3] = dy1
        du[4] = -ω_r1^2*y1+r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m1)
        du[5] = dz2
        du[6] = -ω_z2^2*z2-r_vec[1]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)+V0*(z2/eta0^2*Q2)/m2*sin((ω_zin-ω_zout+delta)*t)
        du[7] = dy2
        du[8] = -ω_r2^2*y2-r_vec[2]*Q1*Q2/(4*pi*epsilon0*(r)^2*m2)
    end
    u0 = [z1_eq+1e-18/sqrt(m1)*vector_z_in[1],0,0,0,z2_eq+1e-18/sqrt(m2)*vector_z_in[2],0,0,0]
    tspan = (0,1000*2*pi/ω_zin);
    problem = ODEProblem(TwoP_COULOMB_DOPPL, u0, tspan)
    solu = solve(problem,Tsit5(),dtmax = 1e-7)
    #in_phase_v = sqrt(137).*(solu[1,:].-z1_eq).*vector_z_in[1]+sqrt(9000/12).*(solu[5,:].-z2_eq).*vector_z_in[2]
    out_phase_v = sqrt(137).*(solu[1,:].-z1_eq).*vector_z_out[1]+sqrt(9000/12).*(solu[5,:].-z2_eq).*vector_z_out[2]
    return minimum(out_phase_v)
end

# optimize(OptimizerFunc,-2*pi*50,2*pi*50,GoldenSection())
plot(p1)