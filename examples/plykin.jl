dt = 5.e-8
function Step(u0::Array{Float64,1},s::Array{Float64,1},n::Int64,
			t0::Float64=0.0)
	u = copy(u0)
	for i = 1:n
		x = u[1]
		y = u[2]
		r2 = x*x + y*y
		t = t0 + (i-1)*dt
		t1 = pi/2.0*t
		st = sin(t1)
		ct = cos(t1)
		cst = cos(pi/4*st)
		sst = sin(pi/4*st)
		cct = cos(pi/4*ct)
		sct = sin(pi/4*ct)
			
		Ω1 = (2.0*x*cct + (1.0 - r2)*sct)/(1 + r2)/(1 + r2)
		Ω2 = (-2.0*x*sst + (1.0 - r2)*cst)/(1 + r2) + 1.0/sqrt(2.0)  
		
		u[1] += dt*(-2.0*s[1]*y*y*Ω1*(cct - x*sct)) + 
		s[2]*y*Ω2*st*(cst - x*sst)
		u[2] += dt*(2*s[1]*y*Ω1*(x*cct + 0.5*(1.0 - x*x + y*y)*sct)) - 
		s[2]*Ω2*st*(x*cst + 0.5*(1.0 - x*x + y*y)*sst)

	end


	return u 
end
