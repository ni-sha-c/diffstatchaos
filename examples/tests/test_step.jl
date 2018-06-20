testName="kuznetsov"
include("../$testName.jl")
using PyPlot

function plot_attractor()

	n_testpoints = 5000
	n_times = 6
	n_stepsi = n_poincare
	u0 = rand(d,n_testpoints)
	u0 .*= (boundaries[d+1:2*d]-boundaries[1:d])
	u0 .+= boundaries[1:d]
	u = zeros(d,n_testpoints,n_times)
	[u[:,:,j] = copy(u0) for j = 1:n_times]
	t = 0.0
	subplot(331)
	plot(u0[1,:],u0[2,:],"o")
	for i = 1:n_testpoints
		u[:,i,1] = Step(u[:,i,1],s0,100*n_stepsi,t)
		t = 100*n_poincare*dt
		for n = 2:n_times
			u[:,i,n] = Step(u[:,i,n-1],s0,100*n_stepsi,t)
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


	r = zeros(n_testpoints)
	theta = zeros(n_testpoints)
	phi = zeros(n_testpoints)
	x = zeros(n_testpoints)
	y = zeros(n_testpoints)
	for i = 1:n_testpoints
		r[i],theta[i],phi[i] = convert_to_spherical(u[:,i,6]) 
		x[i],y[i] = stereographic_projection(u[:,i,6])
	end
	figure()		
	subplot(121)		
	plot(phi,theta,"ko")
	xlim([-pi,pi])
	ylim([0.,pi])
	subplot(122)
	plot(x,y,"ro")


end

function test_tangent()

	n_testpoints = 100
	n_epsi = 8
	
	u0 = rand(d,n_testpoints)
	epsi = logspace(-n_epsi,-1.0,n_epsi)
	vu_fd = zeros(d,n_testpoints,n_epsi)
	vs_fd = zeros(d,n_testpoints,n_epsi)
	vu_ana = zeros(d,n_testpoints)
	vs_ana = zeros(d,n_testpoints)
	u0next = zeros(d)
	v0 = rand(4)
	ds0 = [1.,1.]
	for i = 1:n_testpoints
		u0[:,i] = Step(u0[:,i],s0,10*n_poincare)	
		t = 6.0*rand()
		for k = 1:n_epsi		
			u0next = Step(u0[:,i],s0,1)
			vu_fd[:,i,k] = (Step(u0[:,i] + epsi[k]*v0,
						 s0,1)-u0next)/epsi[k]
			
			vs_fd[:,i,k] = (Step(u0[:,i],s0 + epsi[k]*ds0,1) - 
						u0next)/epsi[k]

		end

		vu_ana[:,i] = tangent_step(v0,u0[:,i],s0,zeros(p))
		vs_ana[:,i] = tangent_step(zeros(d),u0[:,i],s0,ds0)

	
	end		

	erru = zeros(n_epsi)
	errs = zeros(n_epsi)

	[erru[k] = vecnorm(vu_ana[:,:]-vu_fd[:,:,k]) for k = 1:n_epsi]
	[errs[k] = vecnorm(vs_ana[:,:]-vs_fd[:,:,k]) for k = 1:n_epsi]

	figure()
	loglog(epsi,erru)
	figure()
	loglog(epsi,errs)



end

function test_jacobian()


	u0 = rand(4)
	u0[4] *= T
	epsi = 1.e-8
	Jacu = zeros(4,4)
	Jacu[:,1] = (Step(u0 + epsi*[1.0,0.0,0.0,0.0],s0,1) - 
				Step(u0 - epsi*[1.0,0.0,0.0,0.0],s0,1))/
				(2.0*epsi)

  	Jacu[:,2] = (Step(u0 + epsi*[0.0,1.0,0.0,0.0],s0,1) - 
				Step(u0 - epsi*[0.0,1.0,0.0,0.0],s0,1))/
				(2.0*epsi)

	Jacu[:,3] = (Step(u0 + epsi*[0.0,0.0,1.0,0.0],s0,1) - 
				Step(u0 - epsi*[0.0,0.0,1.0,0.0],s0,1))/
				(2.0*epsi)

	Jacu[:,4] = (Step(u0 + epsi*[0.0,0.0,0.0,1.0],s0,1) - 
				Step(u0 - epsi*[0.0,0.0,0.0,1.0],s0,1))/
				(2.0*epsi)

	dFds1 = (Step(u0,s0 + epsi*[1.0,0.0],1)-Step(u0,s0
			- epsi*[1.0,0.0],1))/(2.0*epsi)	


	dFds2 = (Step(u0,s0 + epsi*[0.0,1.0],1)-Step(u0,s0
			- epsi*[0.0,1.0],1))/(2.0*epsi)	


	Jacana = gradFs(u0,s0)
	println(norm(Jacu-Jacana))

	
	v0 = rand(4)
	v0_fd = Jacu*v0 
	v0_hand = tangent_step(v0,u0,s0,zeros(2))
	println(norm(v0_fd - v0_hand))


	v1_fd = v0_fd + dFds1 
	v1_hand = tangent_step(v0,u0,s0,[1.0,0.0])
	println(norm(v1_fd - v1_hand))

	v2_fd = v0_fd + dFds2 
	v2_hand = tangent_step(v0,u0,s0,[0.0,1.0])
	println(norm(v2_fd - v2_hand))


end

#function test_adjoint()

	u0 = rand(4)
	u0[4] *= T
	u1 = Step(u0,s0,1)
	epsi = 1.e-8

	y1 = [0.0,1.,0.,0.]
	y0_ana = adjoint_step(y1,u0,s0,y1)

	y0_fd = zeros(d)
	v0 = zeros(d)
	for i = 1:d
		v0 = zeros(d)
		v0[i] = 1.0
		u0pert = u0 + epsi*v0
		u1pert =  Step(u0pert,s0,1)
		obj2 = u0pert[2] + u1pert[2]

		u0pert = u0 - epsi*v0
		u1pert =  Step(u0pert,s0,1)
		obj1 = u0pert[2] + u1pert[2]

		y0_fd[i] = (obj2 - obj1)/(2.0*epsi) 
	end
	println(norm(y0_fd-y0_ana))

#end

#=function test_testadjoint()
	u = rand(4)
	u[4] *= T
	y1 = rand(4)
	v0 = rand(4)
	v1 = tangent_step(v0,u,s0,zeros(2))
	y0 = adjoint_step(y1,u,s0,zeros(4))
	println(dot(v1,y1))
	println(dot(v0,y0))	
end=#

