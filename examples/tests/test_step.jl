testName="plykin"
include("../$testName.jl")
using PyPlot
s0 = [0.72, 1.9]
#function plot_attractor()

	n_testpoints = 5000
	n_times = 6
	n_steps = 4000
	u0 = -3.0 + 6.0*rand(2,n_testpoints)
	u = zeros(2,n_testpoints,n_times)
	[u[:,:,j] = copy(u0) for j = 1:n_times]
	t = 0.0
	subplot(331)
	plot(u0[1,:],u0[2,:],"o")
	for n = 1:n_times
		for i = 1:n_testpoints
			u[:,i,n] = Step(u[:,i,n],s0,n_steps,t)
		end
		#t += dt*n_steps
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
