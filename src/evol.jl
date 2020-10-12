
"""
	cond_proba_faster(k, mutated_seq, h, J, prob, T = 1)

	"k": position along an amino acid chain
	"mutated_seq": sequence of amino acid of interest
	"h", "J": parameters of DCA Hamiltonian
	"prob": empty vector of length 21
	"T": temperature of DCA

	Returns a vector of length 21 with the conditional probability
	of each aminoacid in position "k" of the sequence given the context.
"""


function cond_proba_faster(k, mutated_seq, h, J, prob, T = 1)
	@inbounds for q_k in 1:21   #k is where the mutation is happening
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
	generate_new_seq_fix_steps_fast(ref_seq, MC_steps, h, J, N, pop, prob, rr, T = 1)

	Official function to use for sampling,
	it is optimized to get "pop", "prob" and "rr"
	in imput from a "generate_msa_xxx" function to spare time.

	Returns a sequence obtained by Gibbs sampling after "MC_steps" steps.
	In input the reference sequence, and all necessary paramters.
"""


function generate_new_seq_fix_steps_fast(ref_seq, MC_steps, h, J, N, pop, prob, rr, T = 1)
	mutated_seq = copy(ref_seq)
	for steps in 1: MC_steps
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = StatsBase.sample!(pop, ProbabilityWeights(cond_proba_faster(pos_mut, mutated_seq, h, J, prob, T)) , rr)[1]
	end
	return mutated_seq
end







# conditional probability

function cond_proba(k, mutated_seq, h, J, T = 1)
	prob = Array{Float64}(undef, 21)
	for q_k = 1:21   #k is where the mutation is happening
		log_proba = h[q_k, k]/T
		for i in 1:202
			if (i < k)
				log_proba += J[mutated_seq[i], q_k , i, k]/T
			else
				log_proba += J[q_k, mutated_seq[i] , k, i]/T
			end
		end
		prob[q_k] = exp(log_proba)
	end
	return normalize(prob,1)
end



"""
	new_amino(pos_mut, mutated_seq, params, flag_p, T = 1)

	Useful only to select for the kind of sampling wanted.
"""


function new_amino(pos_mut, mutated_seq, params, flag_p, T = 1)
	_, _, h_bm, J_bm , fi_tens_pseudo = params
	if flag_p == "p_unif"
		return rand(1:21)
	elseif flag_p == "p_prof"
		return rand(Categorical(fi_tens_pseudo[:, pos_mut]))
	elseif flag_p == "p_cond"
		return rand(Categorical(cond_proba(pos_mut, mutated_seq, h_bm, J_bm, T)))
	else
		error("The name of the probability distribution is inconsistent. Use \"p_unif\", \"p_prof\" or \"p_cond\" ")
	end
end



function generate_msa_fast_fix_steps(output_file, seed_seq, MC_steps, n_seq, h, J, T = 1; n_MSA = 1, prot_name = "TEM-1")
	N = length(seed_seq)
	pop = collect(Int16, 1:21)
	prob = Vector{Float32}(undef, 21)
	rr = Int16.([0])
	@inbounds for i in 1:n_MSA
        file_name = "$output_file.$i"
		FastaWriter(file_name, "w") do file
			@inbounds for j in 1:n_seq
				writeentry(file, "$j |evolved from $prot_name with Gibbs Sampling | $MC_steps MC steps T = $T",
				vec2string(generate_new_seq_fix_steps_fast(seed_seq, MC_steps, h, J, N, pop, prob, rr, T)))
			end
		end
	end
end





# function cond_prob(qk, pos, mutated_seq, h, J)
# 	log_proba = h[q_k, k]
# 	for i in 1:k-1
# 		log_proba += J[mutated_seq[i], q_k ,i, k]
# 	end
# 	for i in k+1:202
# 		log_proba += J[q_k, mutated_seq[i], k, i]
# 	end
# 	prob[q_k] = exp(log_proba/T)
# end
#
#
# function proba(aa_vec, pos, seq, h, J, T)
# 	probab = Array{Float64, 1}(undef, length(aa))
# 	for (i, aa) in enumerate(aa_vec)
# 		probab[i] = cond_prob(aa, seq, h, J, T)
# 	end
# 	return probab
# end




## fix steps



function generate_new_seq_fix_steps_fast_profile(ref_seq, MC_steps, freq)
	mutated_seq = copy(ref_seq)
	for steps in 1: MC_steps
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = rand(Categorical(freq[:, pos_mut]))
	end
	return mutated_seq
end


function generate_new_seq_fix_steps(ref_seq, tot_steps, params, flag_p, T = 1)
	N = length(ref_seq)
	mutated_seq = copy(ref_seq)
	for steps in 1:tot_steps
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = new_amino(pos_mut, mutated_seq, params, flag_p, T)
	end
	return mutated_seq
end



function generate_new_seq_fix_steps_fast2(ref_seq, MC_steps, h, J, flag_p, T = 1, N = 202)
	pop = collect(Int16, 1:21)
	prob = Vector{Float32}(undef, 21)
	mutated_seq = copy(ref_seq)
	muts = 0
	for i in 1:MC_steps
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = StatsBase.sample!(pop, ProbabilityWeights(cond_proba_faster(pos_mut, mutated_seq, h, J, prob, T)) , Int16.([0]))[1]
	end
	return mutated_seq
end

## Fix muts


function generate_new_seq_fast(ref_seq, n_muts, h, J, flag_p, T = 1, N = 202)
	pop = collect(Int16, 1:21)
	prob = Vector{Float32}(undef, 21)
	mutated_seq = copy(ref_seq)
	muts = 0
	while muts < n_muts
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = StatsBase.sample!(pop, ProbabilityWeights(cond_proba_faster(pos_mut, mutated_seq, h, J, prob, T)) , Int16.([0]))[1]
		muts = count_muts(mutated_seq, ref_seq)
	end
	return mutated_seq
end


function generate_new_seq(ref_seq, n_muts, params, flag_p, T=1)
	counts = []
	N = length(ref_seq)
	mutated_seq = copy(ref_seq)
	muts_count = 0
	while muts_count < n_muts
		pos_mut = rand(1:N)
		mutated_seq[pos_mut] = new_amino(pos_mut, mutated_seq, params, flag_p, T)
		muts_count = count_muts(mutated_seq, ref_seq)
        append!(counts, muts_count)
	end
	return mutated_seq, counts
end


function generate_new_seq_fix_muts_fast(ref_seq, n_muts, freq, h, J, flag_p, N, pop, prob, rr; T = 1)
	mutated_seq = copy(ref_seq)
	mm = 0
	while mm < n_muts
		pos_mut = rand(1:N)
		if flag_p == "p_cond"
			mutated_seq[pos_mut] = StatsBase.sample!(pop, ProbabilityWeights(cond_proba_faster(pos_mut, mutated_seq, h, J, prob, T)) , rr)[1]
		elseif flag_p == "p_prof"
			mutated_seq[pos_mut] = rand(Categorical(freq[:, pos_mut]))
		elseif flag_p == "p_rand"
			mutated_seq[pos_mut] = rand(1:21)
		else
			error("No such option for the probability")
		end
		mm = count_muts(mutated_seq, ref_seq)
	end
	return mutated_seq
end


# Generate full msa




function generate_msa(output_file, seed_seq, MC_steps, n_seq, h, J; flag_p = "p_cond", n_MSA = 1)
	for i in 1:n_MSA
		file_name = "$output_file.$i"
		params = (0, 0, h, J, 0)
		FastaWriter(file_name, "w") do file
			for j in 1:n_seq
				seq = generate_new_seq_fix_steps(seed_seq, MC_steps, params, flag_p)
				desc = "$j |evolved from tem1 with $flag_p | $MC_steps MC steps"
				writeentry(file, desc, vec2string(seq))
			end
		end
	end
end



function generate_msa_fix_muts_fast(output_file, seed_seq, muts_vec, params, flag_p; T = 1)
	h = params[3]
	J = params[4]
	freq = params[5]
	N = length(seed_seq)
	pop = collect(Int16, 1:21)
	prob = Vector{Float32}(undef, 21)
	rr = Int16.([0])
	FastaWriter(output_file, "w") do file
		@inbounds for (j, muts) in enumerate(muts_vec)
			writeentry(file, "$j |evolved from tem1 with $flag_p",
			vec2string(generate_new_seq_fix_muts_fast(seed_seq, muts, freq, h, J, flag_p, N, pop, prob, rr)))
		end
	end
end





function generate_msa_fast_persistent(output_path, seed_seq, MC_steps_vec, n_seq, h, J; flag_p = "p_cond")
	N = length(seed_seq)
	pop = Int16.(collect(1:21))
	rr = Int16.([0])
	prob = Vector{Float32}(undef, 21)
	@sync @distributed for j in 1:n_seq
		mutated_seq = copy(seed_seq)
		MC_steps = 1
		while MC_steps <= MC_steps_vec[end]
			pos_mut = rand(1:N)
			mutated_seq[pos_mut] = StatsBase.sample!(pop,ProbabilityWeights(cond_proba_faster(pos_mut, mutated_seq, h, J, prob)), rr)[1]

			if MC_steps in MC_steps_vec
				desc = "$j |evolved from tem1 with $flag_p | $MC_steps MC steps"
				filename = joinpath(output_path, "MC.p_cond.T.1.MC_steps.$MC_steps.n_seqs.$n_seq.pers.fasta" )
				writefasta(filename, [(desc, vec2string(mutated_seq))],  "a")
			end

			MC_steps += 1
		end
	end
end
