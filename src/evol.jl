export evol_msa_fix_steps, stupid_test

stupid_test() = 2


"""
	cond_proba(k, mutated_seq, h, J, prob, T = 1)

	"k": position along an amino acid chain
	"mutated_seq": sequence of amino acid of interest
	"h", "J": parameters of DCA Hamiltonian
	"prob": empty vector of length 21
	"T": temperature of DCA

	Returns a vector of length 21 with the conditional probability of having
	aminoacid " 1 2 3 ... 21" in that sequence context.
"""

function cond_proba(k, mutated_seq, h, J, prob, T = 1)
	@inbounds for q_k in 1:21  
		log_proba = h[q_k, k]
 		for i in 1:k-1
			log_proba += J[mutated_seq[i], q_k ,i, k]
		end
		for i in k+1:202
			log_proba += J[q_k, mutated_seq[i], k, i]
		end
		prob[q_k] = exp(log_proba/T)
	end
	return normalize(prob,1)
end




"""
	evol_seq_fix_steps(ref_seq, MC_steps, h, J, N, pop, prob, rr, T = 1)
	
	
	Returns a sequence obtained by Gibbs sampling after "MC_steps" steps.
	In input the reference sequence, and all necessary paramters.
	Optimized to get "pop", "prob" and "rr" as imput from a "new_msa_xxx" function.

	"ref_seq": amino acid vector of initial sequence
	"MC_steps" : number of Monte Carlo steps performed
	"N" : length of the evolved protein
	"h", "J": parameters of DCA Hamiltonian
	"prob", "prob", "rr": utilities to run code faster
	"T": temperature of DCA
"""

function evol_seq_fix_steps(ref_seq, MC_steps, h, J, N, pop, prob, rr, T = 1)
	mutated_seq = copy(ref_seq)
	@inbounds for steps in 1: MC_steps
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = StatsBase.sample!(pop, ProbabilityWeights(cond_proba(pos_mut, mutated_seq, h, J, prob, T)) , rr)[1]
	end
	return mutated_seq
end





"""
	evol_msa_fix_steps(output_file, seed_seq, MC_steps, n_seq, h, J, T = 1; n_MSA = 1, prot_name = "TEM-1")
		
	Writes a MSA obtained by Gibbs sampling. 

	INPUT:
	"output_file": file where to print MSA in fasta format
	"seed_seq": amino acid vector of initial sequence
	"MC_steps": number of Monte Carlo steps performed to evolve each seq
	"n_seq": number of sequences in the MSA
	"h", "J": parameters of DCA Hamiltonian
	"T": temperature of DCA
	"n_MSA": number of copies of the MSA to be generated
	"prot_name": name of the original protein sequence
"""


function evol_msa_fix_steps(output_file, seed_seq, MC_steps, n_seq; T = 1, n_MSA = 1, prot_name = "TEM-1")
	h, J = extract_par()
	N = length(seed_seq)
	pop = collect(Int16, 1:21)
	prob = Vector{Float32}(undef, 21)
	rr = Int16.([0])
	@inbounds for i in 1:n_MSA
        file_name = "$output_file.$i"
		FastaWriter(file_name, "w") do file
			@inbounds for j in 1:n_seq
				writeentry(file, "$j |evolved from $prot_name with Gibbs Sampling | $MC_steps MC steps, T = $T",
				vec2string(evol_seq_fix_steps(seed_seq, MC_steps, h, J, N, pop, prob, rr, T = 1)))
			end
		end
	end
end






