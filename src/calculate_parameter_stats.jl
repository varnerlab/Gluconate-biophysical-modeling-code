# include -
include("Include.jl")

function compute_some_stats(data_array; dim_index=1)

    mean_val = mean(data_array,dims=dim_index)
    std_val = std(data_array,dims=dim_index)
    return (mean_val,std_val)
end

function main(path_to_ensemble_file::String)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 12.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./src/CellFree.json"

    # load the parameter ensemble file -
    PA = readdlm(path_to_ensemble_file)

    # build the *default* data dictionary -
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    custom_data_dictionary = deepcopy(default_data_dictionary)

    # get elongation rates et al -
    RNAP_elongation_rate = default_data_dictionary["biophysical_constants_dictionary"]["transcription_elongation_rate"]
    Ribsome_elongation_rate = default_data_dictionary["biophysical_constants_dictionary"]["translation_elongation_rate"]
    characteristic_initiation_time_transcription = default_data_dictionary["biophysical_constants_dictionary"]["characteristic_initiation_time_transcription"]
    characteristic_initiation_time_translation = default_data_dictionary["biophysical_constants_dictionary"]["characteristic_initiation_time_translation"]
    default_mRNA_half_life_in_hr = default_data_dictionary["biophysical_constants_dictionary"]["mRNA_half_life_in_hr"]
    default_prot_half_life_in_hr = default_data_dictionary["biophysical_constants_dictionary"]["protein_half_life_in_hr"]

    # get TX_gntr time constant mods -
    tc_mods_tx_ensemble = PA[12,:]
    number_of_genes = 1
    number_of_samples = length(tc_mods_tx_ensemble)

    # compute the time constants for TX -
    tau_tx_ensemble_gntr = zeros(number_of_genes, number_of_samples)
    gene_coding_length_array = custom_data_dictionary["gene_coding_length_array"]

    # what is the length -
    gene_length = gene_coding_length_array[1]

    # what is kE?
    kE = (1/gene_length)*RNAP_elongation_rate
    kI = (1/characteristic_initiation_time_transcription)

    for sample_index = 1:number_of_samples
        tau_factor = (kE/kI)*tc_mods_tx_ensemble[sample_index]
        tau_tx_ensemble_gntr[sample_index] = tau_factor
    end

    # get TX_venus time constant mods -
    tc_mods_tx_ensemble = PA[13,:]
    number_of_genes = 1
    number_of_samples = length(tc_mods_tx_ensemble)

    # compute the time constants for TX -
    tau_tx_ensemble_venus = zeros(number_of_genes, number_of_samples)
    gene_coding_length_array = custom_data_dictionary["gene_coding_length_array"]

    # what is the length -
    gene_length = gene_coding_length_array[2]

    # what is kE?
    kE = (1/gene_length)*RNAP_elongation_rate
    kI = (1/characteristic_initiation_time_transcription)

    for sample_index = 1:number_of_samples
        tau_factor = (kE/kI)*tc_mods_tx_ensemble[sample_index]
        tau_tx_ensemble_venus[sample_index] = tau_factor
    end

    # compute the time constants for TL gntr -
    tc_mods_tl_ensemble = PA[14,:]
    number_of_prots = 1
    number_of_samples = length(tc_mods_tl_ensemble)
    tau_tl_ensemble_gntr = zeros(number_of_prots, number_of_samples)
    prot_coding_length_array = custom_data_dictionary["protein_coding_length_array"]

    # what is the length -
    prot_length = prot_coding_length_array[1]

    # what is kE?
    kE = (1/prot_length)*Ribsome_elongation_rate
    kI = (1/characteristic_initiation_time_translation)

    for sample_index = 1:number_of_samples
        tau_factor = (kE/kI)*tc_mods_tl_ensemble[sample_index]
        tau_tl_ensemble_gntr[sample_index] = tau_factor
    end

    # compute the time constants for TL venus -
    tc_mods_tl_ensemble = PA[15,:]
    number_of_prots = 1
    number_of_samples = length(tc_mods_tl_ensemble)
    tau_tl_ensemble_venus = zeros(number_of_prots, number_of_samples)
    prot_coding_length_array = custom_data_dictionary["protein_coding_length_array"]

    # what is the length -
    prot_length = prot_coding_length_array[2]

    # what is kE?
    kE = (1/prot_length)*Ribsome_elongation_rate
    kI = (1/characteristic_initiation_time_translation)

    for sample_index = 1:number_of_samples
        tau_factor = (kE/kI)*tc_mods_tl_ensemble[sample_index]
        tau_tl_ensemble_venus[sample_index] = tau_factor
    end

    # KL -
    KL_array = Array{Float64,1}()
    for value in vec(PA[22,:])
        push!(KL_array,value)
    end


    # tau_L_1/2 -
    half_life_TL_array = Array{Float64,1}()
    for value in vec(PA[21,:])
        push!(half_life_TL_array, value)
    end

    # dEs -
    dE_array = PA[1:5,:]

    # ns -
    n_array = PA[6:2:10,:]

    # Ks -
    K_array = PA[7:2:11,:]

    # compute half lifes for mRNA_gntr -
    deg_mod_mRNA = PA[16,:]
    number_of_mRNA = 1
    (number_of_samples) = length(deg_mod_mRNA)
    mRNA_half_life_array_gntr = zeros(number_of_samples)
    kD = log(2)/(default_mRNA_half_life_in_hr)
    for sample_index = 1:number_of_samples
        kA = kD*deg_mod_mRNA[sample_index]
        half_life_value = (log(2)/kA)*60    # convert to min
        mRNA_half_life_array_gntr[sample_index] = half_life_value
    end

    # compute half lifes for mRNA_venus -
    deg_mod_mRNA = PA[17,:]
    number_of_mRNA = 1
    (number_of_samples) = length(deg_mod_mRNA)
    mRNA_half_life_array_venus = zeros(number_of_samples)
    kD = log(2)/(default_mRNA_half_life_in_hr)
    for sample_index = 1:number_of_samples
        kA = kD*deg_mod_mRNA[sample_index]
        half_life_value = (log(2)/kA)*60    # convert to min
        mRNA_half_life_array_venus[sample_index] = half_life_value
    end


    # compute half lifes for protein -
    deg_mod_prot = PA[18:20,:]
    (number_of_prot,number_of_samples) = size(deg_mod_prot)
    prot_half_life_array = zeros(number_of_prot, number_of_samples)
    for prot_index = 1:number_of_prot

        kD = log(2)/(default_prot_half_life_in_hr)

        for sample_index = 1:number_of_samples
            kA = kD*deg_mod_prot[prot_index, sample_index]
            half_life_value = (log(2)/kA)*(1/24)    # convert to days
            prot_half_life_array[prot_index, sample_index] = half_life_value
        end
    end

    # gluconate hill equation transfer function detatils
    n_hill = PA[23,:]
    K_hill = PA[24,:]
    
    # return -
    return (n_hill,K_hill,tau_tx_ensemble_gntr,tau_tx_ensemble_venus,tau_tl_ensemble_gntr, tau_tl_ensemble_venus, KL_array, half_life_TL_array, dE_array, mRNA_half_life_array_gntr, mRNA_half_life_array_venus, prot_half_life_array,n_array,K_array)
end

# setup -
path_to_ensemble_file = "simulated/poets_ensemble/PC_T10.dat"

# execute -
(n_hill,K_hill,tau_tx_ensemble_gntr,tau_tx_ensemble_venus, tau_tl_ensemble_gntr,tau_tl_ensemble_venus, KL_array, half_life_TL_array, dE_array, mRNA_half_life_array_gntr, mRNA_half_life_array_venus, prot_half_life_array,n_array,K_array) = main(path_to_ensemble_file)

# compute some stats -
(μ_tx_gntr,σ_tx_gntr) = compute_some_stats(transpose(tau_tx_ensemble_gntr); dim_index = 1)
(μ_tx_venus,σ_tx_venus) = compute_some_stats(transpose(tau_tx_ensemble_venus); dim_index = 1)

(μ_tl_gntr,σ_tl_gntr) = compute_some_stats(transpose(tau_tl_ensemble_gntr); dim_index = 1)
(μ_tl_venus,σ_tl_venus) = compute_some_stats(transpose(tau_tl_ensemble_venus); dim_index = 1)

(μ_KL,σ_KL) = compute_some_stats(KL_array; dim_index=1)
(μ_hl,σ_hl) = compute_some_stats(half_life_TL_array; dim_index=1)

(μ_dE,σ_dE) = compute_some_stats(dE_array; dim_index=2)

(μ_mRNA_half_life_gntr,σ_mRNA_half_life_gntr) = compute_some_stats(mRNA_half_life_array_gntr; dim_index=1)
(μ_mRNA_half_life_venus,σ_mRNA_half_life_venus) = compute_some_stats(mRNA_half_life_array_venus; dim_index=1)


(μ_prot_half_life,σ_prot_half_life) = compute_some_stats(prot_half_life_array; dim_index=2)


(μ_n,σ_n) = compute_some_stats(n_array; dim_index=2)
(μ_K,σ_K) = compute_some_stats(K_array; dim_index=2)

(μ_n_hill,σ_n_hill) = compute_some_stats(n_hill; dim_index=1)
(μ_K_hill,σ_K_hill) = compute_some_stats(K_hill; dim_index=1)


println((μ_K,σ_K))



