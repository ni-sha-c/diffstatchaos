testName="kuznetsov"
include("../$testName.jl")
using PyPlot
#s0 = [0.72, 1.9]
#function plot_attractor()

	n_testpoints = 500
	n_times = 6
	n_steps = n_poincare
	u0 = rand(d,n_testpoints)
	u0 .*= (boundaries[d+1:2*d]-boundaries[1:d])
	u0 .+= boundaries[1:d]
	u = zeros(d,n_testpoints,n_times)
	[u[:,:,j] = copy(u0) for j = 1:n_times]
	t = 0.0
	subplot(331)
	plot(u0[1,:],u0[2,:],"o")
	for i = 1:n_testpoints
		u[:,i,1] = Step(u[:,i,1],s0,100*n_steps,t)
		t = 100*n_poincare*dt
		for n = 2:n_times
			u[:,i,n] = Step(u[:,i,n-1],s0,100*n_steps,t)
			t += 100*n_poincare*dt
		end
		t = 0.0
	end
	subplot(332)
	plot(u[1,:,1],u[2,:,1],"o")
		
	subplot(333)
	plot(u[1,:,2],u[2,:,2],"o")

	subplot(334)
	plot(u[1,:,3],u[2,:,3],"o")

	subplot(335)
	plot(u[1,:,4],u[2,:,4],"o")
		
	subplot(336)
	plot(u[1,:,5],u[2,:,5],"o")



#end
