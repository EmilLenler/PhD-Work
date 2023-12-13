using Plots
using LinearAlgebra
using SparseArrays
using BenchmarkTools


max_n = 1;
T = 0.5*1e-3; #1mK
hbar = 1*1e-34;#hbar
kb = 1.380*1e-23; #Boltzmann constant
omega_i = 272*2*pi*1e3;
omega_o = 772*2*pi*1e3;
beta = 1/(kb*T);
E0 = 1/2*hbar*omega_o;
rho_vector = Vector{Float64}([]);
Z_real = csch(beta*hbar*omega_o/2)/2*csch(beta*hbar*omega_i/2)/2;


for j in range(0,max_n) #loop over in-phase  
    for i in range(0,max_n) #loop over out of phase
    exp_factor = exp(-beta*(hbar*omega_o*(1/2+i)+hbar*omega_i*(1/2+j)));
    push!(rho_vector,exp_factor);
    end #out of phase loop end
end #in phase loop end
print("Z_numeric / Z_analytical: ", sum(rho_vector/Z_real))

rho_vectorNORM = rho_vector/sum(rho_vector)
# println("Numeric Partition / Analytical Parition: ",Z_numeric/Z_real)
rho = Diagonal(rho_vectorNORM);
println(sum(rho_vector))
println(length(rho_vector))
# mean_n = sum(range(0,max_n).*rho_vector/Z_numeric);
# println("Mean n is: ", mean_n)
# println(sizeof(rho))
# println(sizeof(sparse(rho)))
# println(size(rho))
# println(size(range(0,500)))
global sprho = rho;
H = zeros(Complex{Float64},size(rho));
println(sizeof(H))
println(size(rho)[1])
is = [];
os = [];


N = length(range(0,max_n))
for n in range(1,size(rho)[1]) #loop over row
    push!(is,(n-1)÷N)
    push!(os,mod(n-1,N))
    if n < N
        #H[n,n+N+mod(n-1,N)-1] = 1
    elseif n+N > size(rho)[1]
        H[n,n-N+mod(n -1,N)+1] = 1
    else
        #H[n,n+N+mod(n-1,N)-1] = 1
        H[n,n-N+mod(n-1,N)+1] = 1
    end
end
H += transpose(H)
H*=2*pi*500*hbar
println(H)
@btime *(H,sprho)*1im/hbar

rho_0 = sprho;
# println(size(*(H,sprho)-*(sprho,H)))

println(diag((H*sprho-sprho*H)*1im/hbar*1e-6))
# println(diag(*(H,sprho-*(sprho,H))*1im/hbar*1e-6))
# println((diag((*(H,sprho)-*(sprho,H)))*1im/hbar*1e-6)[1])
# @btime for j in 1:1000
#         t = j*1e-3
#         global sprho -= (*(H,sprho)-*(sprho,H))*1im/hbar*1e-3
#         end
#plot(range(1,length(rho_vector)),real(rho_vectorNORM),ls = :dot)