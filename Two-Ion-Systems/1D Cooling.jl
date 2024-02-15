using OrdinaryDiffEq, Plots
using LinearAlgebra
using LaTeXStrings
using DifferentialEquations
using Optim
using Random
using Statistics
hbar = 1 *1e-34

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
Q1 = 1*1.6*1e-19;
omega_z1 = omega_RF*sqrt(-a1/2)
omega_r1 = omega_RF/2*sqrt(q1^2/2+a1)
V_DC = -a1*m1*(3/4*2.7*1e-3)^2*omega_RF^2/(4*Q1*0.248)
V_RF = q1/(2*Q1)*(m1*r0^2*omega_RF^2)
c = 3*1e8
lambda = c/(607.426262*1e12)
alpha = 0.1 #0.5 for radial degeneracy
waist = 500*1e-6
# function Force_PSEUDO(x,y,Q)
#     return [-m1*omega_r1^2*x,-m1*omega_r1^2*y]
# end
function Force_AX(z,DC,Q)::Float64
    return -Q*2*kappa*DC/(z0^2)*z
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
t_end = 20*1e-3
z1_init = sqrt(2*kb*T/(m1*omega_z1^2))
vz1_init = sqrt(2*kb*T/m1)
dt = 2*pi/(omega_RF)/400
println("Initial Time: ", t_init)
println("End Time: ", t_end)
println("Time Step: ", dt)
println("Total number of steps: ",t_end/dt)
println("Initial z1: ", z1_init)
println("Initial z1 velocity: ", vz1_init)
println("Secular z-direction Barium Frequency (kHz):  ", omega_z1/(2*pi*1e3))
state_init = 0
t = t_init
s= 1
z1s = []
ts = []
vz1s = []
y1s = []
vy1s = []
x1s = []
vx1s = []
data_counter = 0
Gamma = 2*pi*15.2*1e6
detuning = -Gamma
k = 2*pi/lambda
kz = k
doppler_detunes = []
I_sat = ()#2*pi^2*Gamma*hbar/(3*lambda^3)
states = []
GBS = []
saturations = linspace(10,0,div(t_end,dt))
while t < t_end
    if t == t_init
        global current_z1 = z1_init
        global current_vz1 = vz1_init
        global state = state_init
    end
    
    ###START OF VERLET BODY
    next_z1 = current_z1 + current_vz1*dt + 1/(2*m1)*Force_AX(current_z1,V_DC,Q1)*dt^2
    next_vz1 = current_vz1 + 1/(2*m1)*(Force_AX(current_z1,V_DC,Q1)+Force_AX(next_z1,V_DC,Q1))*dt
    # RadForce = Force_PSEUDO(current_x1,current_y1,Q1)#Force_RAD(current_x1,current_y1,t,V_DC,V_RF,Q1)
    # next_y1 = current_y1 + current_vy1*dt + 1/(2*m1)*RadForce[2]*dt^2
    # next_x1 = current_x1 + current_vx1*dt + 1/(2*m1)*RadForce[1]*dt^2
    # RadForce_next = Force_PSEUDO(next_x1,next_y1,Q1)#Force_RAD(next_x1,next_y1,t,V_DC,V_RF,Q1)
    # next_vy1 = current_vy1 + 1/(2*m1)*(RadForce[2]+RadForce_next[2])*dt
    # next_vx1 = current_vx1 + 1/(2*m1)*(RadForce[1]+RadForce_next[1])*dt
    ### END OF VERLET BODY
    #Gauss_Beam = (waist/waist_func(current_z1))^2*exp(-2*(current_x1^2+current_y1^2)/(waist_func(current_z1)^2))
    ###START OF DOPPLER/EVENT BODY
    #Gauss_Beam = 1
    if state == 0
        global B12 = s*Gamma*0.5/(1+4*(detuning-(kx*current_vx1+ky*current_vy1+kz*current_vz1))^2/(Gamma^2))
        rnd_nr = rand(Float64)
        if rnd_nr < B12*dt #Absorption
            #next_vx1+=hbar*kx/m1
            #next_vy1+=hbar*ky/m1
            next_vz1+=hbar*kz/m1
            state = 1
        end
    elseif state == 1
        rnd_nr = rand(Float64)
        global B12 = s*0.5*Gamma/(1+4*(detuning-(kx*current_vx1+ky*current_vy1+kz*current_vz1))^2/(Gamma^2))
        if rnd_nr < Gamma*dt #Spont Emission
            #rand_vec = random_unit_vec()
            rand_D = sign(rand(Float64)-0.5)
            #next_vx1 += hbar*k*rand_vec[1]/m1
            #next_vy1 += hbar*k*rand_vec[2]/m1
            next_vz1 += hbar*k*rand_D/m1
            state = 0
        elseif Gamma*dt < rnd_nr
            if rnd_nr < (Gamma+B12)*dt#Stim Emission
                #next_vx1-=kx*hbar/m1
                #next_vy1-=ky*hbar/m1
                next_vz1-=kz*hbar/m1
                state = 0
            end
        end
    end
    ###END OF DOPPLER/EVENT BODY
    
    if mod(data_counter,500) == 0
        push!(z1s,current_z1)
        push!(vz1s,current_vz1)
        push!(ts,t)
        # push!(y1s,current_y1)
        # push!(vy1s,current_vy1)
        # push!(x1s,current_x1)
        # push!(vx1s,current_vx1)
        data_counter = 0
        #push!(GBS,Gauss_Beam)
        #push!(doppler_detunes,(detuning+(kx*current_vx1+ky*current_vy1+kz*current_vz1))/Gamma)
        #push!(states,state)
        #println(B12*dt)
    end
    current_z1 = next_z1
    current_vz1 = next_vz1
    # current_y1 = next_y1
    # current_vy1 = next_vy1
    # current_x1 = next_x1
    # current_vx1 = next_vx1
    global t += dt
    global data_counter +=1
end
vs = LinRange(-1000*sqrt(2*5*kb*t/m1),1000sqrt(2*5*kb*t/m1),1000)
E_kin = 1/2*m1.*(vz1s.^2)
U_pot = 1/2*m1.*(z1s.^2 .*omega_z1^2)
p1 = plot(ts*1e3,z1s*1e6,linewidth = 2,xlabel = "Time (ms)",ylabel = "\$z_1\$ (µm)")
#p7 = plot(ts,GBS)
p8 = plot(ts*1e3,(E_kin .+ U_pot)./(kb/2) * 1e3,xlabel = "Time (ms)", ylabel = "Kinetic Energy/1.5kb (mK)",ylimits = (0,5),hline = 1.5)
p9 = plot(ts*1e3,vz1s,xlabel = "Time (ms)", ylabel = "z-velocity")
hline!(p8,[0.5])
savefig(p1,"z_plot.png")
savefig(p8,"total_energy")
savefig(p9,"z velocity")
