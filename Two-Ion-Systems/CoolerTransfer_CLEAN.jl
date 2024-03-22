using OrdinaryDiffEq, Plots
using LinearAlgebra
using LaTeXStrings
using DifferentialEquations
using Optim
using Random
using Statistics
hbar = 1 *1e-34

function print_section()
    println("---------------------------------------------------------------------------------------------------------------------")
end
Random.seed!(95757) #Seed for reproducibility
T = 10 * 1e-3
kb = 1.38*1e-23
kappa = 0.248;
z0 = 3/4*2.7*1e-3;
r0 = 3/4*3.5*1e-3;
a1 = -0.001;
q1 = 0.21801801801801804;
omega_RF = 2*pi*5.2*1e6;
m1 = 135*1.66*(10^(-27))
m2 = 9000*1.66*(10^(-27))
Q1 = 1*1.6*1e-19;
Q2 = 0#24*1.6*1e-19;
omega_z1 = omega_RF*sqrt(-a1/2)
omega_r1 = omega_RF/2*sqrt(q1^2/2+a1)
V_DC = -a1*m1*(3/4*2.7*1e-3)^2*omega_RF^2/(4*Q1*0.248)
V_RF = q1/(2*Q1)*(m1*r0^2*omega_RF^2)
c = 3*1e8
lambda = c/(607.426262*1e12)
alpha = 0.1 #0.5 for radial degeneracy
waist = 500*1e-6
epsilon0 = 8.85*1e-12

rho = Q2/Q1
mu = m2/m1

a2 = rho/mu * a1
q2 = rho/mu * q1

omega_z2 = omega_RF*sqrt(-a2/2)

omega_pond_1 = Q1*V_RF/(sqrt(2)*omega_RF*m1*r0^2)
omega_x1 = sqrt(omega_pond_1^2 - alpha*omega_z1^2)
omega_y1 = sqrt(omega_pond_1^2-(1-alpha)*omega_z1^2)

omega_pond2 = Q2*V_RF/(sqrt(2)*omega_RF*m2*r0^2)
omega_x2 = sqrt(omega_pond2^2-alpha*omega_z2^2)
omega_y2 = sqrt(omega_pond2^2-(1-alpha)*omega_z2^2)

zeq1 = cbrt(Q1*Q2/(4*pi*epsilon0*m1*omega_z1^2)*(1/(1+1/rho)^2))
zeq2 = -1/rho * zeq1
V_tickle = 0#0.5*V_DC
print_section()
println("DC Voltage: ",V_DC)
println("RF Voltage: ", V_RF )
println("Tickle Voltage: ", V_tickle)
print_section()
println("Equilibrium positions (µm): " , "z1 = ", zeq1 * 1e6, "    z2 = " ,zeq2 * 1e6)


deltaZ = zeq1-zeq2
K11 = omega_z1^2 + Q1*Q2/(4*pi*epsilon0*m1) * (2/deltaZ^3)
K12 = -Q1*Q2/(4*pi*epsilon0*sqrt(m1*m2)) * (2/deltaZ^3)
K22 = omega_z2^2 + Q1*Q2/(4*pi*epsilon0*m2) * (2/deltaZ^3)

K33 = omega_x1^2 - Q1*Q2/(4*pi*epsilon0*m1) * (1/deltaZ^3)
K34 = -0.5*K12
K44 = omega_x2^2 - Q1*Q2/(4*pi*epsilon0*m2) * (1/deltaZ^3)

K55 = omega_y1^2 - Q1*Q2/(4*pi*epsilon0*m1) * (1/deltaZ^3)
K56 = K34
K66 = omega_y2^2 - Q1*Q2/(4*pi*epsilon0*m2) * (1/deltaZ^3)


KZ = [K11 K12
     K12 K22]

KX = [K33 K34
      K34 K44]

KY = [K55 K56
      K56 K66]

z_freqs = sqrt.(eigvals(KZ))
x_freqs = sqrt.(eigvals(KX))
y_freqs = sqrt.(eigvals(KY))

z_vecs = eigvecs(KZ)
x_vecs = eigvecs(KX)
y_vecs = eigvecs(KY)
print_section()
println("z frequencies are (kHz) : ", (z_freqs[1])/(2*pi*1e3),"     ", (z_freqs[2])/(2*pi*1e3))
println("z eigenvectors are: ", z_vecs[1,:], "        ", z_vecs[2,:])
print_section()
println("x frequencies are (kHz) : ", (x_freqs[1])/(2*pi*1e3),"     ", (x_freqs[2])/(2*pi*1e3))
println("x eigenvectors are: ", x_vecs[1,:], "        ", x_vecs[2,:])
print_section()
println("y frequencies are (kHz) : ", (y_freqs[1])/(2*pi*1e3),"     ", (y_freqs[2])/(2*pi*1e3))
println("y eigenvectors are: ", y_vecs[1,:], "        ", y_vecs[2,:])
print_section()
function Force_RAD(x,y,t,DC,RF,Q)
    return [2*alpha*Q*kappa*DC/(z0^2)*x + Q*RF/(r0^2)*cos(omega_RF*t)*x,2*Q*kappa*DC/(z0^2)*(1-alpha)*y-Q*RF/(r0^2)*cos(omega_RF*t)*y]
end
function Force_PSEUDO1(x,y)
    return [-m1*omega_x1^2*x,-m1*omega_y1^2*y]
end
function Force_PSEUDO2(x,y)
    return [-m2*omega_x2^2*x,-m2*omega_y2^2*y]
end
function Force_AX(z,DC,Q)::Float64
    return -Q*2*kappa*DC/(z0^2)*z
end
function Force_COULOMB(x1,x2,y1,y2,z1,z2,chrg1,chrg2)
    r = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)
    return chrg1*chrg2/(4*pi*epsilon0).*[(x1-x2)/(r^3),(y1-y2)/r^3,(z1-z2)/r^3]
end
function random_unit_vec()
    u = rand(Float64)
    v = rand(Float64)
    phi = 2*pi*u
    theta = acos(1-2*v)
    v_x = sin(theta)*cos(phi)
    v_y = sin(theta)*sin(phi)
    v_z = cos(theta)
    return [v_x,v_y,v_z]
end
t_init = 0
t_end = 1*1e-3
z1_init = zeq1#+1e-18/sqrt(m1)*z_vecs[1,1]#-sqrt(2*kb*T/(m1*omega_z1^2))
r1_init = 0#sqrt(2*kb*T/(m1*omega_r1^2))
vz1_init = 0#sqrt(2*kb*T/m1)
y1_init = 0#-r1_init
vy1_init = 0#sqrt(2*kb*T/m1)
x1_init = 0#r1_init
vx1_init = 0#-sqrt(2*kb*T/m1)

z2_init = zeq2#+1e-18/sqrt(m2)*z_vecs[1,2]#sqrt(2*kb*T/(m2*omega_z2^2))
y2_init = 0#sqrt(2*kb*T/(m2*omega_y2^2))
x2_init = 0#sqrt(2*kb*T/(m2*omega_x2^2))
vz2_init = 0#-sqrt(2*kb*T/m2)
vx2_init = 0#vz2_init
vy2_init = 0#-vz2_init


dt = 2*pi/(omega_RF)/400*(pi/3)
println("Initial Time: ", t_init)
println("End Time: ", t_end)
println("Time Step: ", dt)
println("Total number of steps: ",floor(t_end/dt))
print_section()
println("Initial z1: ", z1_init)
println("Initial z1 velocity: ", vz1_init)
println("Initial x1: ", x1_init)
println("Initial x1 velocity: ", vx1_init)
println("Initial y1: ", y1_init)
println("Initial y1 velocity: ", vy1_init)
print_section()
println("Initial z2: ", z2_init)
println("Initial z2 velocity: ", vz2_init)
println("Initial x2: ", x2_init)
println("Initial x2 velocity: ", vx2_init)
println("Initial y2: ", y2_init)
println("Initial y2 velocity: ", vy2_init)
print_section()
println("Secular z-direction Barium Frequency (kHz):  ", omega_z1/(2*pi*1e3))
println("Secular x-direction Barium Frequency (kHz):  ", omega_y1/(2*pi*1e3))
println("Secular y-direction Barium Frequency (kHz):  ", omega_x1/(2*pi*1e3))
print_section()
state_init = 0
t = t_init
s = 0.1
z1s = []
ts = []
vz1s = []
y1s = []
vy1s = []
x1s = []
vx1s = []

z2s = []
vz2s = []
y2s = []
vy2s = []
x2s = []
vx2s = []

proj_zi = []
proj_zo = []
if sign(z_vecs[1,1]) == sign(z_vecs[1,2])
    z_i = z_vecs[1,:]
    if sign(z_i[1]) == -1
        z_i .*= -1
    end
    z_o = z_vecs[2,:]
else
    z_i = z_vecs[2,:]
    if sign(z_i[1]) == -1
        z_i .*= -1
    end
    z_o = z_vecs[1,:]
end



data_counter = 1
Gamma = 2*pi*15.2*1e6
detuning = -Gamma/2
k = 2*pi/lambda
kx = 1/2*k
ky = 1/2*k
kz = 1/sqrt(2)*k
doppler_detunes = []
I_sat = ()#2*pi^2*Gamma*hbar/(3*lambda^3)
states = []
GBS = []
while t < t_end
    if t == t_init
        global current_z1 = z1_init
        global current_vz1 = vz1_init
        global current_y1 = y1_init
        global current_vy1 = vy1_init
        global current_x1 = x1_init
        global current_vx1 = vx1_init

        global current_z2 = z2_init
        global current_vz2 = vz2_init
        global current_y2 = y2_init
        global current_vy2 = vy2_init
        global current_x2 = x2_init
        global current_vx2 = vx2_init

        global state = state_init
    end
    
    ###START OF VERLET BODY
    Fc_now = Force_COULOMB(current_x1,current_x2,current_y1,current_y2,current_z1,current_z2,Q1,Q2)
    RadForce1 = Force_PSEUDO1(current_x1,current_y1)
    RadForce2 = Force_PSEUDO2(current_x2,current_y2)
    next_z1 = current_z1 + current_vz1*dt + 1/(2*m1)*(Force_AX(current_z1,V_DC,Q1))#+Force_AX(current_z1,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q1)+Fc_now[3])*dt^2
    next_y1 = current_y1 + current_vy1*dt + 1/(2*m1)*(RadForce1[2]+Fc_now[2])*dt^2
    next_x1 = current_x1 + current_vx1*dt + 1/(2*m1)*(RadForce1[1]+Fc_now[1])*dt^2
    next_z2 = current_z2 + current_vz2*dt + 1/(2*m2)*(Force_AX(current_z2,V_DC,Q2))#+Force_AX(current_z2,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q2)-Fc_now[3])*dt^2
    next_y2 = current_y2 + current_vy2*dt + 1/(2*m2)*(RadForce2[2]-Fc_now[2])*dt^2
    next_x2 = current_x1 + current_vx2*dt + 1/(2*m2)*(RadForce2[1]-Fc_now[1])*dt^2

    Fc_next = Force_COULOMB(next_x1,next_x2,next_y1,next_y2,next_z1,next_z2,Q1,Q2)
    RadForce1_next = Force_PSEUDO1(next_x1,next_y1)
    RadForce2_next = Force_PSEUDO2(next_x2,next_y2)

    next_vz1 = current_vz1 + 1/(2*m1)*(Force_AX(current_z1,V_DC,Q1)+Force_AX(next_z1,V_DC,Q1)+Force_AX(current_z1,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q1)+Force_AX(next_z1,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q1)+Fc_now[3]+Fc_next[3])*dt
    next_vy1 = current_vy1 + 1/(2*m1)*(RadForce1[2]+RadForce1_next[2]+Fc_now[2]+Fc_next[2])*dt
    next_vx1 = current_vx1 + 1/(2*m1)*(RadForce1[1]+RadForce1_next[1]+Fc_now[1]+Fc_next[1])*dt
    next_vz2 = current_vz2 + 1/(2*m2)*(Force_AX(current_z2,V_DC,Q2)+Force_AX(next_z2,V_DC,Q2)+Force_AX(current_z2,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q2)+Force_AX(next_z2,V_tickle*cos((z_freqs[1]-z_freqs[2])*t),Q2)-Fc_now[3]-Fc_next[3])*dt
    next_vy2 = current_vy2 + 1/(2*m2)*(RadForce2[2]+RadForce2_next[2]-Fc_now[2]-Fc_next[2])*dt
    next_vx2 = current_vx2 + 1/(2*m2)*(RadForce2[1]+RadForce2_next[1]-Fc_now[1]-Fc_next[1])*dt
    ### END OF VERLET BODY
    #Gauss_Beam = (waist/waist_func(current_z1))^2*exp(-2*(current_x1^2+current_y1^2)/(waist_func(current_z1)^2))
    ###START OF DOPPLER/EVENT BODY
    #Gauss_Beam = 1
    if state == 0
        global B12 = s*Gamma*0.5/(1+4*(detuning-(kx*current_vx1+ky*current_vy1+kz*current_vz1))^2/(Gamma^2))
        rnd_nr = rand(Float64)
        if rnd_nr < B12*dt #Absorption
            next_vx1+=hbar*kx/m1
            next_vy1+=hbar*ky/m1
            next_vz1+=hbar*kz/m1
            state = 1
        end
    elseif state == 1
        rnd_nr = rand(Float64)
        global B12 = s*0.5*Gamma/(1+4*(detuning-(kx*current_vx1+ky*current_vy1+kz*current_vz1))^2/(Gamma^2))
        if rnd_nr < Gamma*dt #Spont Emission
            rand_vec = random_unit_vec()
            next_vx1 += hbar*k*rand_vec[1]/m1
            next_vy1 += hbar*k*rand_vec[2]/m1
            next_vz1 += hbar*k*rand_vec[3]/m1
            state = 0
        elseif rnd_nr < (Gamma+B12)*dt
            #if rnd_nr < (Gamma+B12)*dt#Stim Emission
            next_vx1-=kx*hbar/m1
            next_vy1-=ky*hbar/m1
            next_vz1-=kz*hbar/m1
            state = 0
        end
    end
    ###END OF DOPPLER/EVENT BODY
    
    if mod(data_counter,500) == 0
        push!(z1s,current_z1)
        push!(vz1s,current_vz1)
        push!(ts,t)
        push!(y1s,current_y1)
        push!(vy1s,current_vy1)
        push!(x1s,current_x1)
        push!(vx1s,current_vx1)
        push!(z2s,current_z2)
        push!(vz2s,current_vz2)
        push!(y2s,current_y2)
        push!(vy2s,current_vy2)
        push!(x2s,current_x2)
        push!(vx2s,current_vx2)
        push!(proj_zi,(current_z1-zeq1)*z_i[1]*sqrt(135)+(current_z2-zeq2)*z_i[2]*sqrt(9000))
        push!(proj_zo,(current_z1-zeq1)*z_o[1]*sqrt(135)+(current_z2-zeq2)*z_o[2]*sqrt(9000))
    end
    current_z1 = next_z1
    current_vz1 = next_vz1
    current_y1 = next_y1
    current_vy1 = next_vy1
    current_x1 = next_x1
    current_vx1 = next_vx1

    current_z2 = next_z2
    current_vz2 = next_vz2
    current_y2 = next_y2
    currenty_vy2 = next_vy2
    current_x1 = next_x1
    current_vx1 = next_vx1
    global t += dt
    global data_counter +=1
end
E_kin = 1/2*m1.*(vz1s.^2 .+ vx1s.^2 .+vy1s.^2)
U_pot = 1/2*m1.*(z1s.^2 .*omega_z1^2 .+ omega_r1.^2 .*(y1s.^2 .+ x1s.^2))
p1 = plot(ts*1e3,z1s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$z_1\$ (µm)",label = "\$z_1\$")
p2 = plot(ts*1e3,y1s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$y_1\$ (µm)",label = "\$y_1\$")
p3 = plot(ts*1e3,x1s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$x_1\$ (µm)",label = "\$x_1\$")
p4 = plot(ts .*1e3, 1/2*m1*omega_z1^2/kb .*z1s.^2 .*1e3,xlabel = "Time (ms)",ylabel = "Potential energy in z (mK)",ylimits = (0,2))
p5 = plot(ts .*1e3, 1/2*m1*omega_r1^2/kb .*y1s.^2 .*1e3,xlabel = "Time (ms)",ylabel = "Potential energy in y (mK)",ylimits = (0,2))
p6 = plot(ts .*1e3, 1/2*m1*omega_r1^2/kb .*x1s.^2 .*1e3,xlabel = "Time (ms)",ylabel = "Potential energy in x (mK)",ylimits = (0,2))
p9 = plot(ts*1e3,vz1s,xlabel = "Time (ms)", ylabel = L"\dot{z}_1" * "  (m/s)",label = L"\dot{z}_1")
p10 = plot(ts*1e3,vy1s,xlabel = "Time (ms)", ylabel = L"\dot{y}_1" * "  (m/s)",label = L"\dot{y}_1")
p11 = plot(ts*1e3,vx1s,xlabel = "Time (ms)", ylabel = L"\dot{x}_1" * "  (m/s)", label = L"\dot{x}_1")
savefig(p1,"z1_plot.png")
savefig(p2,"y1_plot.png")
savefig(p3,"x1_plot.png")
savefig(p4,"z1_pot_in_mK")
savefig(p5,"y1_pot_in_mK")
savefig(p6,"x1_pot_in_mK")
savefig(p9,"z1 velocity")
savefig(p10,"y1 velocity")
savefig(p11,"x1 velocity")
println(z_vecs[1,1],z_vecs[1,2])

z2pl = plot(ts*1e3,z2s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$z_2\$ (µm)")
x2pl = plot(ts*1e3,x2s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$x_2\$ (µm)")
y2pl = plot(ts*1e3,y2s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$y_2\$ (µm)")
savefig(z2pl,"z2_plot")
savefig(x2pl,"x2_plot")
savefig(y2pl,"y2_plot")


zpl_scale = maximum(vcat(proj_zi,proj_zo))
println(zpl_scale)
zipl = plot(ts*1e3,proj_zi,ylimits = (-2*zpl_scale,2*zpl_scale),label = "In-phase motion")
zopl = plot(ts*1e3,proj_zo,xlabel = "Time (ms)",ylimits = (-2*zpl_scale,2*zpl_scale), label ="Out of phase motion")
zpl = plot(zipl,zopl, layout = (2,1))
savefig(zpl, "Both Z Projections")

#Do a mean E_kin calculation:
L = length(E_kin)
mean_Es = []
interval_lgth = 10
mean_ts = []
for j in 1:Integer(floor(L/interval_lgth))
    E_k_mean = sum(E_kin[1+(j-1)*interval_lgth : 1+j*interval_lgth])/interval_lgth
    push!(mean_Es,E_k_mean)
    push!(mean_ts,sum(ts[1+(j-1)*interval_lgth: 1+j*interval_lgth])/interval_lgth)
end
t_avg = ts[interval_lgth+1]-ts[1]
t_avg_str = string(t_avg*1e6)
p12 = plot(mean_ts*1e3,mean_Es/(1.5*kb)*1e3,ylimits = (0,2),xlabel = "Time (ms)", ylabel = "Energy / 1.5kb (mK)",title = "Averaging time: "*t_avg_str* " (µs)")
savefig(p12,"Temperature from velocity mean (mK)")