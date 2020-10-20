"""

"""
function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)
    N = size(h)[2]
    q = size(h)[1]
    E = 0.0
    for i = 1:N
        E -= h[S[i],i]
        for j = (i+1):N
			E -= J[S[i],S[j],i,j]
		end
	end
return E
end






# takes the parameters of francesco and changes in a way that we get a tensor (for J)
# and that the J[a, b, i, j] is such that  1 == "A", ... "21" == "-"
function read_par_BM(path::AbstractString, q, N)
	data = readdlm(path,' ', use_mmap = true)[:, 2:6]
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = data[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = data[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end


