import Distributions
import StatsBase
import JSON

ITER = 50000
BURNIN = 1000
THIN = 25
BETA_PROPOSAL_SHAPE = 55
INFILE = ARGS[1]
OUTFILE = ARGS[2]

function mean(x)
    return sum(x)/length(x)
end

function step_x(x)
    d = Distributions.Beta(max(1, BETA_PROPOSAL_SHAPE*x), max(1, BETA_PROPOSAL_SHAPE*(1-x)))
    return rand(d)
end

function l_dbinom(k, n, p)
#    println([k, n, p])
    d = Distributions.Binomial(n,p)
    return Distributions.logpdf(d, k)
end

function move_pos_in_fams_bound_0size(x, fam_size)
    xstep = x + StatsBase.sample([-1, 0, 1])
    if xstep < 0
        return 0
    elseif xstep > fam_size
        return fam_size
    else
        return xstep
    end
end

function get_f1_ct_from_parent_ct(parent_ct, fam_size)
    if parent_ct == 0
        return 0
    elseif parent_ct == 4
        return fam_size
    else
        d = Distributions.BetaBinomial(fam_size, parent_ct, 4-parent_ct)
        return rand(d)
    end
end

function get_prob(mu_val, f2_rr_ct, parent_pos_ct, f1_pos_ct, f1_ct, f2_ct)
#this is the core bit that computes the likelihood
#parent and f1 is tracking alleles (2N pop size) while f2 is tracking N pop size, b/c trait is recessive
#assumes mu is scalar and rest of variable are vectors
#so when they toher terms are scalars, like "x", need to convert to "[x]" to make it work
    lnlk = 0
    for idx in 1:length(f2_rr_ct)
        parent_freq = parent_pos_ct[idx] / 4. #parent R alelle freq
        f1_freq = f1_pos_ct[idx]/f1_ct[idx]   #f1 R allele freq
        lnlk += l_dbinom(parent_pos_ct[idx], 4., mu_val) # parent gen likelihood
        lnlk += l_dbinom(f1_pos_ct[idx], f1_ct[idx], parent_freq) # f1 gen likelihood
        lnlk += l_dbinom(f2_rr_ct[idx], f2_ct[idx], f1_freq*f1_freq) #f2 likelihood, when trait recessive E[survivors] is f1_freq^2, for domimnant it'd be 1-f1_freq^2
    end
    return lnlk
end

function update_theta(curr_prob_l, prop_prob_l, curr_theta, prop_theta)
#generic func. to do acceptance sampling on log prob called theta
#returns new theta and 0/1 for reject/accept tracking
    if prop_prob_l > curr_prob_l
        return (prop_theta, 1) #acceptance
    else
        l_prob_move = prop_prob_l - curr_prob_l
        l_rand = log(rand())
        if l_rand < l_prob_move
            return(prop_theta, 1) #acceptance
        else
            return(curr_theta, 0) #reject proposal
        end
    end
end

xfile = []
open(INFILE) do file
    for ln in eachline(file)
        z = split(ln, "\t")
        push!(xfile, z[1:3])
    end
end

f1_ct, f2_ct, rr_ct = ([],[],[])
for idx in 2:length(xfile) #exclude header
    m = map(u->parse(Int, u), xfile[idx,])
    append!(f1_ct, m[1])
    append!(f2_ct, sum(m[2:3]))
    append!(rr_ct, m[3])
end

curr_mu = 0.1 #starting val, this seems fine

curr_parent_pos_ct = zeros(length(f1_ct))
curr_f1_pos_ct = zeros(length(f1_ct))

for idx in 1:length(rr_ct)
    if rr_ct[idx] > 0
        curr_parent_pos_ct[idx] = 1
        curr_f1_pos_ct[idx] = ceil(f1_ct[idx]*0.25)
    end
end

pos = length([i for i in rr_ct if i > 0])
neg = length(rr_ct)*4 - pos
println([pos, neg])

mu_posterior, parent_pos_posterior, f1_pos_posterior, acceptance_ct, acceptance_par_and_f1 = ([], [], [], [], [])

for i in 1:ITER
    if i >= BURNIN
        if i % THIN == 0 #past burnin and thin, so record
            append!(mu_posterior, curr_mu)
            push!(parent_pos_posterior, curr_parent_pos_ct)
            push!(f1_pos_posterior, curr_f1_pos_ct)
            if length(acceptance_ct) > 0
                mean_acc = mean(acceptance_ct)
                println(ITER-i, " ", round(mean(mu_posterior), digits=3), " ", round(mean_acc, digits=3), " ", length(mu_posterior))
                println(map((xa)-> round(xa, digits=3), StatsBase.quantile(mu_posterior, [1/40, 1/2, 39/40])), "\n")
            end
        end
    end
    
    prop_mu = step_x(curr_mu)
    curr_prob_l = get_prob(curr_mu, rr_ct, curr_parent_pos_ct, curr_f1_pos_ct, f1_ct, f2_ct)
    prop_prob_l = get_prob(prop_mu, rr_ct, curr_parent_pos_ct, curr_f1_pos_ct, f1_ct, f2_ct)
    
    global curr_mu, accepted = update_theta(curr_prob_l, prop_prob_l, curr_mu, prop_mu)
    if i >= BURNIN
        append!(acceptance_ct, accepted) #tracking acceptance rate
    end
    
    accpt = []
    for idx in 1:length(f1_ct)
        #parent and f1 correlated, so move together
        prop_parent_pos_ct = copy(curr_parent_pos_ct)
        prop_parent_pos_ct[idx] = move_pos_in_fams_bound_0size(curr_parent_pos_ct[idx], 4)
        prop_f1_pos_ct = copy(curr_f1_pos_ct)
        prop_f1_pos_ct[idx] = get_f1_ct_from_parent_ct(prop_parent_pos_ct[idx], f1_ct[idx])
        
        curr_prob_l = get_prob(curr_mu, [rr_ct[idx]], [curr_parent_pos_ct[idx]], [curr_f1_pos_ct[idx]], [f1_ct[idx]], [f2_ct[idx]])
        prop_prob_l = get_prob(curr_mu, [rr_ct[idx]], [prop_parent_pos_ct[idx]], [prop_f1_pos_ct[idx]], [f1_ct[idx]], [f2_ct[idx]])
        
        curr_parent_and_f1_vec, accp_val = update_theta(curr_prob_l, prop_prob_l, [curr_parent_pos_ct, curr_f1_pos_ct], [prop_parent_pos_ct, prop_f1_pos_ct])
        global curr_parent_pos_ct = curr_parent_and_f1_vec[1]
        global curr_f1_pos_ct = curr_parent_and_f1_vec[2]
        push!(accpt, accp_val)
    end
    push!(acceptance_par_and_f1, accpt)
end

D = Dict("mu"=>mu_posterior, "parent_pos"=> parent_pos_posterior, "f1_pos"=>f1_pos_posterior, "f1_ct"=>f1_ct, "rr_ct"=>rr_ct)

open(OUTFILE, "w") do f
    JSON.print(f, D)
end
println(StatsBase.quantile(mu_posterior, [1/40, 1/4, 1/2, 3/4, 39/40]))
