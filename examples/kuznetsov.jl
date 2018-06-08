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
function Step(u0::Array{Float64,1},s::Array{Float64,1},n::Int64,
			 t0::Float64)

	u = copy(u0)

	

	for i = 1:n

		x = u[1]
		y = u[2]
		z = u[3]
		σ = 0.0
		a = 0.0	
		r2 = x^2 + y^2 + z^2	
		r = sqrt(r2)
		
		t = (t0 + (i-1)*dt)%T
		σ = diff_rot_freq(t)
		a = rot_freq(t)		
	
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

function tangent_step(v0::Array{Float64,1},u::Array{Float64,1},
					 s::Array{Float64,1},
					 ds::Array{Float64,1},
					 t::Float64)

	x = u[1]
	y = u[2]
	z = u[3]

	dx = v0[1]
	dy = v0[2]
	dz = v0[3]
	v = copy(v0)

	σ = 0.0
	a = 0.0	
	r2 = x^2 + y^2 + z^2	
	r = sqrt(r2)
		
	t = t%T
	σ = diff_rot_freq(t)
	a = rot_freq(t)

	
	coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
	coeff2 = s[1]*(1. - σ*σ - a*a)
	coeff3 = s[2]*a*a*(1.0 - r)		

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
				coeff3*dx)  

	v[2] += dt*(coeff1*0.5*dx + 
				dcoeff1_dz*0.5*x*dz + 
				dcoeff2_ds1*y*x*x*ds[1] + 
				coeff2*dy*x*x + 
				coeff2*2.0*x*dx*y + 
				dcoeff3_ds2*ds[2]*y + 
				dcoeff3_dx*y*dx + 
				dcoeff3_dy*y*dy +
		 		dcoeff3_dz*y*dz +
			 	coeff3*dy) 

	v[3] += dt*(-0.5*a*pi*dx + 
				dcoeff3_ds2*z*ds[2] + 
				dcoeff3_dx*z*dx + 
				dcoeff3_dy*z*dy + 
				dcoeff3_dz*z*dz + 
				coeff3*dz)	 

	return v



end

function rot_freq(t::Float64)

	if(t >= 2.0 && t < 3.0)
		return -1.0
	end

	if(t >= 5.0 && t < 6.0)
		return 1.0
	end
	
	return 0.0

end

function diff_rot_freq(t::Float64)

	if(t >= 1.0 && t < 2.0)
		return -1.0
	end

	if(t >= 4.0 && t < 5.0)
		return 1.0
	end

	return 0.0


end

function adjoint_step(y1::Array{Float64,1},u::Array{Float64,1},
					 s::Array{Float64,1},
					 dJ::Array{Float64,1},t::Float64)


	y0 = copy(y1)

	x = u[1]
	y = u[2]
	z = u[3]
	#t = u[4]
	
	r2 = x^2 + y^2
	r = sqrt(r2)
	σ = diff_rot_freq(t)
	a = rot_freq(t)		
	coeff1 = σ*pi*0.5*(z*sqrt(2) + 1)
	coeff2 = s[1]*(1. - σ*σ - a*a)
	coeff3 = s[2]*a*a*(1.0 - r)		

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
		

	return y0

end	
