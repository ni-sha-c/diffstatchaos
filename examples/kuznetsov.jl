dt = 1.e-2
s0 = [1.0,1.0]
d = 3
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
			
		if(t >= 1.0 && t < 2.0)
			σ = -1.0
		elseif(t >= 4.0 && t < 5.0)
			σ = 1.0
		elseif(t >= 2.0 && t < 3.0)
			a = -1.0
		elseif(t >= 5.0 && t < 6.0)
			a = 1.0
		end	
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


