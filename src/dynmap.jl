#module dynmap

include("solenoid_contract.jl")
using PGFPlots
s = [1.4, 0.0]
u0 = [1.,1.,1.]
d = 3

function plot_trajectories()

	n_samples = 1
	n_steps = 20000
	ut = zeros(d,n_steps,n_samples)
	p = plot()
	for k = 1:n_samples
		ut[:,1,k] = rand(d)
		for n = 1:n_steps-1
			ut[:,n+1,k] = Step(ut[:,n,k],s,1)
		end 
		for n = 5000:100:n_steps-1 
			#plot!(p,ut[1,n,k],ut[2,n,k],
			#ut[3,n,k],
			#linestyle=:dot,
			#linewidth=1,label="",
			#xlabel="x",ylabel="y",zlabel="z")
			#sleep()
		end
	end
	#default(show = true)	
	#return p
end

function plot_objective()
  
	n_samples = 100
	n_steps = 10000
	sk = copy(s)
	Φ = zeros(n_samples)
	for k = 1:n_samples
		sk[1] = sk[1] + (k-1)/n_samples*2.0 
		u = rand(d)
		u = Step(u, sk, 5000)
		
		Φ[k] = 0
		for n = 1:n_steps
			
			u = Step(u, sk, 1)
			Φ[k] += objective(u,sk)/n_steps
		end
	end

	p = plot(linspace(1.,1.+2.*(n_samples-1)/n_samples,n_samples),
		 Φ,legend=false,
		 annotate=(2.5,80,text("Objective",:blue)),xlabel="R",linewidth=3,
		 linestyle=:dash,linecolor=:blue)

	return p
end  

function verify24()

	n_more_samples = 1000
	expΦ0 = 0.
	expΦoF0 = 0.
	rhs2 = 0.0
	rhs3 = 0.0
	expΦs = 0.
	expΦoFs = 0.
	ss = copy(s)
	ss[1] += 0.005
	n_samples = 1000
	uarr0 = rand(d, n_samples)
	uarrs = rand(d, n_samples)
	for k = 1:n_samples
		uarr0[:,k] = Step(uarr0[:,k], s, 10000)
		uarrs[:,k] = Step(uarrs[:,k], ss, 10000)
	end
	for k = 1:n_more_samples
		println("k= ", k)	
		for i = 1:n_samples
			uarr0[:,i] = Step(uarr0[:,i], s, 100)
			expΦ0 += objective(uarr0[:,i],s)/n_samples
			rhs2 += compute242(uarr0[:,i],s)/n_samples
			ΦoF = objective(Step(uarr0[:,i],s,1),s)
			expΦoF0 +=  ΦoF/n_samples 
		
		end
		for i = 1:n_samples
			uarrs[:,i] = Step(uarrs[:,i], ss, 100)	
			expΦs += objective(uarrs[:,i],ss)/n_samples
			uarrsTmp = Step(uarrs[:,i], s, 1)
			expΦoFs +=  objective(uarrsTmp,ss)/n_samples 
		end
	
	end	
	
	lhs = (expΦs-expΦ0)/0.005/n_more_samples
	println(lhs)
	rhs1 = (expΦoFs-expΦoF0)/0.005/n_more_samples
	rhs2 /= n_more_samples
	println(rhs1+rhs2)
	

end

function Ψs()

	n_samples = 1000
	n_steps = 5000
	uarr = rand(d,n_samples)
	for k = 1:n_samples
		uarr[:,k] = Step(uarr[:,k],s,n_steps)
	end 
	Ψsv0 = 0.
	ΨsMsv0 = 0.
	expv0oFsdiv∇Fsinv = 0.	
	for k = 1:n_samples

		println("k = ", k)
		Ψsv0 += -1.0*Divv(uarr[:,k],s)/n_samples
		ΨsMsv0 += -1.0*compute29RHS1(uarr[:,k],s)/n_samples
		expv0oFsdiv∇Fsinv += compute29RHS2(uarr[:,k],s)/n_samples


	end

	lhs = Ψsv0
	rhs = ΨsMsv0 + expv0oFsdiv∇Fsinv

	println(lhs)
	println(rhs)
	

end
function decompose(v,u,w0)

	rdir = [u[1],u[2],0.0]
	zdir = [0.0, 0.0, 1.0]
	w1 = zdir
	w2 = rdir - dot(rdir,w1)*w1
	w2 /= norm(w2)
	vu = dot(v,w0)*w0/(norm(w0)^2.0)		
	vs = v - vu
					
	
	return vs,vu
end



function stableSens()

	n_iter = 1000
	n_samples = 1000
	n_steps = 5000
	uarr = rand(d,n_samples)

	for k = 1:n_samples
		uarr[:,k] = Step(uarr[:,k],s,n_steps)
	end 
	v = zeros(3)
	v0 = rand(3)
	dΦds_stable = zeros(n_samples)
	for i = 1:n_samples

		println("i= ",i)
		v = zeros(3)
		u = uarr[:,i]
		v0 = rand(3)
		for j = 1:n_iter
			dFds = ∂F∂s(u,s)
			nablaFs = gradFs(u,s)
			
			u = Step(u,s,1)
			gradφ = nablaφ(u,s)    
			v0 = nablaFs*v0
			v0 /= norm(v0)
			
			v = nablaFs*v + dFds
			v, dFds_unstable = decompose(v, u, v0)
			dΦds_stable[i] += dot(gradφ, v)/n_iter
		end
	

	end
	a = Axis(Plots.Histogram(dΦds_stable, bins=20), ymin=0, 
    xlabel=L"$\frac{d\Phi}{ds}$", title="Sensitivity under
	stable perturbations")

	save("stable_sensitivity.pdf",a)
	
end
#function compute_sensitivity()

	s0 = [1.4,0.0]
	s = copy(s0)
	p = 2
	K = 10
	Kprime = 10
	n_samples = 5000000
	n_runup = 1000

	n_steps = n_samples + K +
	Kprime + 1 
			
	
	n_bins_t = 50
	n_bins_r = 50
	rmax = 3.0
	binsize_r = rmax/n_bins_r
	binsize_t = 2*pi/n_bins_t
	u = rand(d)
	
	u = Step(u,s,n_runup)	
	b = zeros(p)'
	dΦds_stable = zeros(n_bins_t,n_bins_r,p)
	dΦds_unstable = zeros(n_bins_t,n_bins_r,p)
	
	φ = zeros(n_bins_t,n_bins_r)
	
	g = zeros(K,p)
	ϕhat = zeros(d)'
	v0 = rand(d)
	v = zeros(d,p)
	
	t = linspace(-pi,pi-binsize_t,n_bins_t)
	r = linspace(0., rmax - binsize_r, n_bins_r)
		
	dFds = zeros(d,p)
	dFds_unstable = zeros(d,p)
	nablaFsinv = zeros(d,d)
	gsum_mean = zeros(p)
	gsum_history = zeros(n_steps,p)
	for n = 1:n_steps

		dFds = ∂F∂s(u,s) 
		nablaFs = gradFs(u,s)
		b = ∇Fsinvcolon∂Fs∂s(u,s)
		q = div∇Fsinv(u,s)
		ϕhat = -1.0*q + 
			transpose(\(nablaFs', ϕhat'))
		v0 = nablaFs*v0
		v0 /= norm(v0)
		ϕhat = dot(ϕhat', v0) * v0'
		u = Step(u,s,1)
		for ip=1:p
			
			g[(n-1)%K+1,ip] = ϕhat*dFds_unstable[:,ip] - b[ip]   
			gsum_mean[ip] += sum(g[:,ip])
			v[:,ip] = nablaFs*v[:,ip] + dFds[:,ip]
			gsum_history[n,ip] = sum(g[:,ip])
			v[:,ip], dFds_unstable[:,ip] = decompose(v[:,ip], u, v0)
			if(n>=K+Kprime+1)
			
				for binno_t=1:n_bins_t
					for binno_r=1:n_bins_r
						φ[binno_t,binno_r] = objective(u,s,t[binno_t],binsize_t,
											r[binno_r],binsize_r)	
						dΦds_unstable[binno_t,binno_r,ip] += 
								#φ[binno_t,binno_r]*(sum(g[:,ip])
								#		-gsum_mean[ip]/n)/n_samples
								φ[binno_t,binno_r]*(sum(g[:,ip]))/n_samples

						gradφ = nablaφ(u,s,t[binno_t],binsize_t,
								r[binno_r],binsize_r)
						dΦds_stable[binno_t,binno_r,ip] += dot(gradφ, v[:,ip])/n_samples
					end	
				end
			end
			
		end
	end
	
#=
#	ds = 0.05
#	u1 = rand(3,p)
#	u = rand(3)
#	u = Step(u,s0,5000)

#	n_fd_samples = n_samples*500
	ρs = zeros(n_bins_t,n_bins_r,p)
	ρ0 = zeros(n_bins_t,n_bins_r)


		
#	for ip = 1:p
#		s = copy(s0)
#		s[ip] = s0[ip] + ds
#		u1[:,ip] = Step(u1[:,ip],s,5000)
#	end
	
#	φ = zeros(n_bins)
	


	for n = 1:n_fd_samples
		println(n)
 		for ip = 1:p
			s = copy(s0)	
			s[ip] = s0[ip] + ds	
			u = Step(u,s0,1)
			u1[:,ip] = Step(u1[:,ip],s,1)
		
			for binno_t = 1:n_bins_t
				for binno_r = 1:n_bins_r
					ρ0[binno_t,binno_r] += objective(u,
							s0,t[binno_t],r[binno_r])/n_fd_samples
					ρs[binno_t,binno_r,ip] += objective(u1[:,ip],s,
						t[binno_t],r[binno_r])/n_fd_samples 
				end			
			end

		end
	end
=#	


	t1 = repmat(t,1,n_bins_r)
	t1 = t1[:]
	r1 = repmat(r',n_bins_t,1)
	r1 = r1[:]
	dΦds_1 = dΦds_unstable[:,:,1] + dΦds_stable[:,:,1]
	writedlm("variables/50/dPhids_r_alg.txt",dΦds_1);
	drhot_r_alg = sum(dΦds_1,2); 
	dΦds_1 = dΦds_1[:]
	dΦds_2 = dΦds_unstable[:,:,2] + dΦds_stable[:,:,2]
	writedlm("variables/50/dPhids_t_alg.txt",dΦds_2);
	drhot_t_alg = sum(dΦds_2,2); 
	dΦds_2 = dΦds_2[:]
	sens_plot = Axis(Plots.Linear3(t1,r1,dΦds_1),style="only marks",
				xlabel=L"$\theta$",ylabel="r")
	push!(sens_plot,Plots.Linear3(t1,r1,dΦds_2))

	push!(sens_plot_rhot_t, Plots.Linear(t + zeros(n_bins_t),drhot_t_alg[:]))
	#save("drhot_t_comp.pdf",sens_plot_rhot_t)
	
	
	
	#dρds_1 = (ρs[:,:,1]-ρ0[:,:])/ds
	#dρds_1 = dρds_1[:]
	#dρds_2 = (ρs[:,:,2]-ρ0[:,:])/ds
	#dρds_2 = dρds_2[:]
	#push!(sens_plot,Plots.Linear3(t,r,dρds_1))	
	#push!(sens_plot,Plots.Linear3(t,r,dρds_2))


	#save("dSRBds.pdf",sens_plot)

#end




#end
