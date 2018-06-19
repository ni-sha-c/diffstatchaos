dt = 1.e-2
s0 = [1.0,1.0]
d = 3
p = 2
boundaries = ones(2*d)
boundaries[1] = -1.0
boundaries[2] = -1.0
boundaries[3] = -1.0
T = 6.0
n_poincare = ceil(Int64,T/dt)
function Step(u0::Array{Float64,1},s::Array{Float64,1},n::Int64)

	u = copy(u0)

	

	for i = 1:n

		x = u[1]
		y = u[2]
		z = u[3]
		σ = 0.0
		a = 0.0	
		r2 = x^2 + y^2 + z^2	
		r = sqrt(r2)
		
		σ = diff_rot_freq(u[4])
		a = rot_freq(u[4])		
	
		coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
		coeff2 = s[1]*(1. - σ*σ - a*a)
		coeff3 = s[2]*a*a*(1.0 - r)		

		u[1] += dt*(-1.0*coeff1*y - 
				coeff2*x*y*y + 
				0.5*a*pi*z + coeff3*x)

		u[2] += dt*(coeff1*0.5*x + 
				coeff2*y*x*x + 
				coeff3*y)

		u[3] += dt*(-0.5*a*pi*x + coeff3*z)

		u[4] = (u[4] + dt)%T

	end 


	return u
end

function convert_to_spherical(u::Array{Float64,1})

	x = u[1]	
	y = u[2]	
	z = u[3]
	r = sqrt(x^2 + y^2 + z^2)
	θ = acos(z/r)
	ϕ = atan2(y,x)	
	return r,θ,ϕ

end

function stereographic_projection(u::Array{Float64,1})

	x = u[1]
	y = u[2]
	z = u[3]
	deno = x + z + sqrt(2.0)

	re_part = (x - z)/deno
	im_part = y*sqrt(2.0)/deno

	return re_part,im_part

end
function tangent_source(v0::Array{Float64,1}, u::Array{Float64,1},
						s::Array{Float64,1}, ds::Array{Float64,1})


	v = copy(v0)
	v[1] += dt*(-1.0*dcoeff2_ds1*ds[1]*x*y*y + 
				dcoeff3_ds2*ds[2]*x)
	v[2] += dt*(dcoeff2_ds1*ds[1]*y*x*x + 
				dcoeff3_ds2*ds[2]*y)
	v[3] += dt*dcoeff3_ds2*ds[2]*z
	

	return v

end
function ∂F∂s(u::Array{Float64,1},s::Array{Float64,1})

	dfds = zeros(d,p)
	ds1 = [1.0, 0.0]
	ds2 = [0.0, 1.0]
	dfds[:,1] = tangent_source(zeros(d),u,s,ds1)
	dfds[:,2] = tangent_source(zeros(d),u,s,ds2)
	return dfds

end

function gradFs(u::Array{Float64,1},s::Array{Float64,1})

	x = u[1]
	y = u[2]
	z = u[3]
	t = u[4]

	σ = 0.0
	a = 0.0	
	r2 = x^2 + y^2 + z^2	
	r = sqrt(r2)
		
	t = t%T
	
	σ = diff_rot_freq(t)
	a = rot_freq(t)
	dσ_dt = ddiff_rot_freq_dt(t)
	da_dt = drot_freq_dt(t)

	coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
	coeff2 = s[1]*(1. - σ*σ - a*a)
	coeff3 = s[2]*a*a*(1.0 - r)		

	dcoeff1_dt = pi*0.5*(z*sqrt(2) + 1)*dσ_dt
	dcoeff2_dt = s[1]*(-2.0)*(σ*dσ_dt + a*da_dt)
	dcoeff3_dt = s[2]*(1.0 - r)*2.0*a*da_dt


	dcoeff1_dz = σ*pi*0.5*sqrt(2)
	dcoeff2_ds1 = coeff2/s[1]
	dcoeff3_ds2 = coeff3/s[2]
	dcoeff3_dx = s[2]*a*a*(-x)/r
	dcoeff3_dy = s[2]*a*a*(-y)/r
	dcoeff3_dz = s[2]*a*a*(-z)/r
		
	dFds = zeros(d,d)

	dFdu[1,1] = 1.0 + dt*(-coeff2*y*y +
					coeff3 + dcoeff3_dx)

	dFdu[1,2] = 
	
	
	
	-1.0*dcoeff1_dz*y*dz - 1.0*
				coeff1*dy - dcoeff2_ds1*ds[1]*x*y*y - 
			 - coeff2*x*2.0*y*dy + 
				0.5*a*pi*dz + dcoeff3_ds2*ds[2]*x + 
				dcoeff3_dx*x*dx + 
				dcoeff3_dy*x*dy + 
				dcoeff3_dz*x*dz +
				coeff3*dx - 1.0*dcoeff1_dt*y*dtime - 
				dcoeff2_dt*x*y*y*dtime + 
				0.5*da_dt*pi*z*dtime + 
				dcoeff3_dt*x*dtime)

	v[2] += dt*(coeff1*0.5*dx + 
				dcoeff1_dz*0.5*x*dz + 
				dcoeff2_ds1*y*x*x*ds[1] + 
				coeff2*dy*x*x + 
				coeff2*2.0*x*dx*y + 
				dcoeff3_ds2*ds[2]*y + 
				dcoeff3_dx*y*dx + 
				dcoeff3_dy*y*dy +
		 		dcoeff3_dz*y*dz +
			 	coeff3*dy + dcoeff1_dt*0.5*x*dtime -
				dcoeff2_dt*y*x*x*dtime + 
				dcoeff3_dt*y*dtime) 

	v[3] += dt*(-0.5*a*pi*dx + 
				dcoeff3_ds2*z*ds[2] + 
				dcoeff3_dx*z*dx + 
				dcoeff3_dy*z*dy + 
				dcoeff3_dz*z*dz + 
				coeff3*dz - 
				0.5*pi*x*da_dt*dtime + 
				dcoeff3_dt*z*dtime)
	 

	return v



end


function tangent_step(v0::Array{Float64,1},u::Array{Float64,1},
					 s::Array{Float64,1},
					 ds::Array{Float64,1})

	x = u[1]
	y = u[2]
	z = u[3]
	t = u[4]
	dx = v0[1]
	dy = v0[2]
	dz = v0[3]
	dtime = v0[4]
	v = copy(v0)

	σ = 0.0
	a = 0.0	
	r2 = x^2 + y^2 + z^2	
	r = sqrt(r2)
		
	t = t%T
	
	σ = diff_rot_freq(t)
	a = rot_freq(t)
	dσ_dt = ddiff_rot_freq_dt(t)
	da_dt = drot_freq_dt(t)


	
	coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
	coeff2 = s[1]*(1. - σ*σ - a*a)
	coeff3 = s[2]*a*a*(1.0 - r)		

	dcoeff1_dt = pi*0.5*(z*sqrt(2) + 1)*dσ_dt
	dcoeff2_dt = s[1]*(-2.0)*(σ*dσ_dt + a*da_dt)
	dcoeff3_dt = s[2]*(1.0 - r)*2.0*a*da_dt


	dcoeff1_dz = σ*pi*0.5*sqrt(2)
	dcoeff2_ds1 = coeff2/s[1]
	dcoeff3_ds2 = coeff3/s[2]
	dcoeff3_dx = s[2]*a*a*(-x)/r
	dcoeff3_dy = s[2]*a*a*(-y)/r
	dcoeff3_dz = s[2]*a*a*(-z)/r
		


	v[1] += dt*(-1.0*dcoeff1_dz*y*dz - 1.0*
				coeff1*dy - dcoeff2_ds1*ds[1]*x*y*y - 
				coeff2*y*y*dx - coeff2*x*2.0*y*dy + 
				0.5*a*pi*dz + dcoeff3_ds2*ds[2]*x + 
				dcoeff3_dx*x*dx + 
				dcoeff3_dy*x*dy + 
				dcoeff3_dz*x*dz +
				coeff3*dx - 1.0*dcoeff1_dt*y*dtime - 
				dcoeff2_dt*x*y*y*dtime + 
				0.5*da_dt*pi*z*dtime + 
				dcoeff3_dt*x*dtime)

	v[2] += dt*(coeff1*0.5*dx + 
				dcoeff1_dz*0.5*x*dz + 
				dcoeff2_ds1*y*x*x*ds[1] + 
				coeff2*dy*x*x + 
				coeff2*2.0*x*dx*y + 
				dcoeff3_ds2*ds[2]*y + 
				dcoeff3_dx*y*dx + 
				dcoeff3_dy*y*dy +
		 		dcoeff3_dz*y*dz +
			 	coeff3*dy + dcoeff1_dt*0.5*x*dtime -
				dcoeff2_dt*y*x*x*dtime + 
				dcoeff3_dt*y*dtime) 

	v[3] += dt*(-0.5*a*pi*dx + 
				dcoeff3_ds2*z*ds[2] + 
				dcoeff3_dx*z*dx + 
				dcoeff3_dy*z*dy + 
				dcoeff3_dz*z*dz + 
				coeff3*dz - 
				0.5*pi*x*da_dt*dtime + 
				dcoeff3_dt*z*dtime)
	 

	return v



end

function rot_freq(t::Float64)
	
	a0 = -1.0
	a1 = 0.0
	a2 = 1.0

	c0 = 2.0
	c1 = 3.0
	c2 = 5.0
	c3 = 6.0 
	c4 = 0.0

	slope = 20.0
	est = exp(slope*t)
	esc0 = exp(slope*c0)
	esc1 = exp(slope*c1)
	esc2 = exp(slope*c2)
	esc3 = exp(slope*c3)
	esc4 = exp(slope*c4)

	fn0 = (a1*esc0 + a0*est)/(esc0 + est)	
	fn1 = (a0*esc1 + a1*est)/(esc1 + est)
	fn2 = (a1*esc2 + a2*est)/(esc2 + est)
	fn3 = (a2*esc3 + a1*est)/(esc3 + est)
	fn4 = (a2*esc4 + a1*est)/(esc4 + est)

	return fn0 + fn1 + fn2 + fn3 + fn4




end

function diff_rot_freq(t::Float64)

	a0 = -1.0
	a1 = 0.0
	a2 = 1.0

	c0 = 1.0
	c1 = 2.0
	c2 = 4.0
	c3 = 5.0 

	slope = 20.0
	est = exp(slope*t)
	esc0 = exp(slope*c0)
	esc1 = exp(slope*c1)
	esc2 = exp(slope*c2)
	esc3 = exp(slope*c3)
	
	fn0 = (a1*esc0 + a0*est)/(esc0 + est)	
	fn1 = (a0*esc1 + a1*est)/(esc1 + est)
	fn2 = (a1*esc2 + a2*est)/(esc2 + est)
	fn3 = (a2*esc3 + a1*est)/(esc3 + est)

	return fn0 + fn1 + fn2 + fn3

end

function ddiff_rot_freq_dt(t::Float64)
	
	a0 = -1.0
	a1 = 0.0
	a2 = 1.0

	c0 = 1.0
	c1 = 2.0
	c2 = 4.0
	c3 = 5.0 

	slope = 20.0
	est = exp(slope*t)
	esc0 = exp(slope*c0)
	esc1 = exp(slope*c1)
	esc2 = exp(slope*c2)
	esc3 = exp(slope*c3)


	dfn0 = esc0*est*slope*(a0-a1)/(esc0 + est)/(esc0 + est)
	dfn1 = esc1*est*slope*(a1-a0)/(esc1 + est)/(esc1 + est)
	dfn2 = esc2*est*slope*(a2-a1)/(esc2 + est)/(esc2 + est)
	dfn3 = esc3*est*slope*(a1-a2)/(esc3 + est)/(esc3 + est)

	return dfn0 + dfn1 + dfn2 + dfn3 


end

function drot_freq_dt(t::Float64)
	
	a0 = -1.0
	a1 = 0.0
	a2 = 1.0

	c0 = 2.0
	c1 = 3.0
	c2 = 5.0
	c3 = 6.0 
	c4 = 0.0

	slope = 20.0
	est = exp(slope*t)
	esc0 = exp(slope*c0)
	esc1 = exp(slope*c1)
	esc2 = exp(slope*c2)
	esc3 = exp(slope*c3)
	esc4 = exp(slope*c4)

	fn0 = (a1*esc0 + a0*est)/(esc0 + est)	
	fn1 = (a0*esc1 + a1*est)/(esc1 + est)
	fn2 = (a1*esc2 + a2*est)/(esc2 + est)
	fn3 = (a2*esc3 + a1*est)/(esc3 + est)
	fn4 = (a2*esc4 + a1*est)/(esc4 + est)

	dfn0 = esc0*est*slope*(a0-a1)/(esc0 + est)/(esc0 + est)
	dfn1 = esc1*est*slope*(a1-a0)/(esc1 + est)/(esc1 + est)
	dfn2 = esc2*est*slope*(a2-a1)/(esc2 + est)/(esc2 + est)
	dfn3 = esc3*est*slope*(a1-a2)/(esc3 + est)/(esc3 + est)
	dfn4 = esc4*est*slope*(a1-a2)/(esc4 + est)/(esc4 + est)

	return dfn0 + dfn1 + dfn2 + dfn3 + dfn4


end


function adjoint_step(y1::Array{Float64,1},u::Array{Float64,1},
					 s::Array{Float64,1},
					 dJ::Array{Float64,1})


	y0 = copy(y1)

	x = u[1]
	y = u[2]
	z = u[3]
	t = u[4]
	
	r2 = x^2 + y^2
	r = sqrt(r2)
	σ = diff_rot_freq(t)
	a = rot_freq(t)	
	dσdt = ddiff_rot_freq_dt(t)
	dadt = drot_freq_dt(t)
	coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
	coeff2 = s[1]*(1. - σ*σ - a*a)
	coeff3 = s[2]*a*a*(1.0 - r)		

	dcoeff1dt = pi*0.5*(z*sqrt(2) + 1)*dσdt
	dcoeff2dt = s[1]*(-2.0)*(σ*dσdt + a*dadt)
	dcoeff3dt = s[2]*(1.0 - r)*2.0*a*dadt

	dcoeff1dz = σ*pi*0.5*sqrt(2)
	dcoeff3dx = s[2]*a*a*(-x)/r 
	dcoeff3dy = s[2]*a*a*(-y)/r
	dcoeff3dz = s[2]*a*a*(-z)/r 

	y0[1] += y1[1]*dt*(-1.0*coeff2*y*y) + 
			y1[1]*dt*coeff3 + 
			y1[1]*dt*x*dcoeff3dx + 
			y1[2]*dt*coeff1*0.5 + 
			y1[2]*dt*y*2.0*x + 
			y1[2]*dt*dcoeff3dx*y + 
			y1[3]*dt*(-0.5)*a*pi + 
			y1[3]*dt*z*dcoeff3dx 	

	y0[2] += y1[1]*dt*(-1.0)*coeff1 - 
			y1[1]*dt*(-1.0)*coeff2*x*2.0*y + 
			y1[1]*dt*dcoeff3dy*x + 
			y1[2]*dt*coeff2*x*x + 
			y1[2]*dt*coeff3 + 
			y1[2]*dt*dcoeff3dy*y + 
			y1[3]*dt*dcoeff3dy*z

	
	y0[3] += y1[1]*dt*(-1.0)*dcoeff1dz*y + 
			y1[1]*dt*0.5*a*pi + 
			y1[1]*dt*dcoeff3dz*x + 
			y1[2]*dt*dcoeff1dz*0.5*x + 
			y1[2]*dt*y*dcoeff3dz + 
			y1[3]*dt*z*dcoeff3dz + 
			y1[3]*dt*coeff3


	y0[4] += -1.0*y1[1]*dt*dcoeff1dt*y - 
			 y1[1]*dt*x*y*y*dcoeff2dt + 
			 y1[1]*dt*x*dcoeff3dt +
			 y1[1]*dt*0.5*pi*z*dadt +  
			 y1[2]*dt*0.5*x*dcoeff1dt + 
			 y1[2]*dt*y*x*x*dcoeff2dt +
			 y1[2]*dt*y*dcoeff3dt + 
			 y1[3]*dt*dcoeff3dt*z + 
			 y1[3]*dt*(-0.5)*pi*x*dadt 
			 

		

	return y0

end	
