include("Include.jl")
# mean center -
function mean_center_array(results_array::Array{Float64,2})::Array{Float64,2}

    # get the size -
    (NR,NC) = size(results_array)
    scaled_array = zeros(NR,NC)

    for col_index = 1:NC

        data_col = results_array[:,col_index]
        mu_value = mean(data_col)
        std_value = std(data_col)

        for row_index = 1:NR
            scaled_array[row_index,col_index] = (data_col[row_index] - mu_value)/(std_value)
        end
    end

    return scaled_array
end

# computes the model performance -
function model_performance(parameter_guess_array,index)

    # what is the host_type?
    host_type = :cell_free

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 12.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./src/CellFree.json"
    #path_to_data_dir = "$(pwd())/data"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = deepcopy(default_data_dictionary)

   # Phase 1: parameter update =========================================================================== #
   # update the paramaters in the model data dictionary -
   # for now - lets only search over dG's -
   R = model_data_dictionary["R"]
   T_K = model_data_dictionary["T_K"]

   # compute W -
   tmp_W_array = Float64[]
   for index = 1:5
       parameter_guess = parameter_guess_array[index]
       value = exp(-1*(parameter_guess))
       push!(tmp_W_array,value)
   end

   # update the control W's -
   control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
   control_parameter_dictionary["W_GntR_RNAP"] = tmp_W_array[1]
   control_parameter_dictionary["W_GntR_sigma_70"] = tmp_W_array[2]
   control_parameter_dictionary["W_Venus_RNAP"] = tmp_W_array[3]
   control_parameter_dictionary["W_Venus_sigma_70"] = tmp_W_array[4]
   control_parameter_dictionary["W_Venus_GntR"] = tmp_W_array[5]
   model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
   #print(control_parameter_dictionary)

   binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
   binding_parameter_dictionary["n_GntR_sigma_70"] = parameter_guess_array[6]
   binding_parameter_dictionary["K_GntR_sigma_70"] = parameter_guess_array[7]
   binding_parameter_dictionary["n_Venus_sigma_70"] = parameter_guess_array[8]
   binding_parameter_dictionary["K_Venus_sigma_70"] = parameter_guess_array[9]
   binding_parameter_dictionary["n_Venus_GntR"] = parameter_guess_array[10]
   binding_parameter_dictionary["K_Venus_GntR"] = parameter_guess_array[11]
   model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
   #print(binding_parameter_dictionary)
   # time constant modifier -

   time_constant_modifier_array = [
       0.0							;	# 1	GntR
       0.0							;	# 2	Venus
       0.0							;	# 3	sigma_70
       parameter_guess_array[12]	;	# 4	mRNA_GntR
       parameter_guess_array[13]	;	# 5	mRNA_Venus
       1.0							;	# 6	mRNA_sigma_70
       parameter_guess_array[14]	;	# 7	protein_GntR
       parameter_guess_array[15]	;	# 8	protein_Venus
       1.0							;	# 9	protein_sigma_70
   ]

   model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

   # setup degradation_modifier_array -

   degradation_modifier_array = [
       0.0	;	# 1	GntR
       0.0	;	# 2	Venus
       0.0	;	# 3	sigma_70
       parameter_guess_array[16]	;	# 4	mRNA_GntR
       parameter_guess_array[17]	;	# 5	mRNA_Venus
       1.0	;	# 6	mRNA_sigma_70
       parameter_guess_array[18]	;	# 7	protein_GntR
       parameter_guess_array[19]	;	# 8	protein_Venus
       parameter_guess_array[20]	;	# 9	protein_sigma_70
   ]

   model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array

   # update the translation time -
   model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[21]

   # lastly, update KL -
   biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
   biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[22]
   model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

   # gluconate GntR binding parameters
   gluconate_parameter_dictionary = model_data_dictionary["gluconate_parameter_dictionary"]
   gluconate_parameter_dictionary["n_gluconate_GntR"] = parameter_guess_array[23]
   gluconate_parameter_dictionary["K_gluconate_GntR"] = parameter_guess_array[24]
   model_data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary

   # update the transcription capacity parameters
   model_data_dictionary["transcription_capacity_delay"] = parameter_guess_array[25]
   model_data_dictionary["transcription_capacity_slope"] = parameter_guess_array[26]
   
   # update the translation capacity parameters
   model_data_dictionary["translation_capacity_delay"] = parameter_guess_array[27]
   model_data_dictionary["translation_capacity_slope"] = parameter_guess_array[28]


   # grab defaults -
   species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
   protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
   gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
   time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
   initial_condition_array = model_data_dictionary["initial_condition_array"]

   # # get gene IC -
   idx_gene = findall(x->x==:gene,species_symbol_type_array)
   gene_abundance_array = initial_condition_array[idx_gene]

   # Precompute the translation parameters -
   translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
   model_data_dictionary["translation_parameter_array"] = translation_parameter_array

   # Precompute the kinetic limit of transcription -
   transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
   model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array


   # Dilution degrdation matrix -
   dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
   model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

   # ===================================================================================================== #
   #print(model_data_dictionary)
   # Phase 2:  solve model equations ===================================================================== #
   # solve the balance equations -
   (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
    # ===================================================================================================== #

    # Phase 3: compute the model performance metrics ====================================================== #
    p_Venus_AUC = integrate(TSIM,XSIM[:,index])
    # ===================================================================================================== #

    # return the performance_array -
    return p_Venus_AUC
end

function main(path_to_ensemble_file::String,index)

    # setup the sensitivity function -
    SF(P) = model_performance(P,index)

    # setup ranges -
    sample_bounds_array = Array{Tuple,1}()
    ensemble_array = readdlm(path_to_ensemble_file)
    (number_of_parameters,number_of_trials) = size(ensemble_array)
    for parameter_index = 1:(number_of_parameters)

        # get row of parameters -
        parameter_row = ensemble_array[parameter_index,:]
        min_value = minimum(parameter_row)
        max_value = maximum(parameter_row)

        # create the tuple -
        tmp_tuple = (min_value,max_value)

        # cache -
        push!(sample_bounds_array,tmp_tuple)
    end

    @show index

    # do the global sensitivity analysis -
    sensitivity_results = GlobalSensitivity.gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),sample_bounds_array)

    # return -
    return sensitivity_results
end

# setup paths -
path_to_ensemble_file = "./src/simulated/poets_ensemble/PC_T10.dat"

# compute a sensitivity array for the AUC of each species -
species_index_array = [5 8 4 7] 
number_of_species = length(species_index_array)
number_of_parameters = 28
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array

    global results_array

    # conduct senstivity analysis -
    sensitivity_results = main(path_to_ensemble_file,species_index)

    # get the μ and σ^2
    mu = sensitivity_results.means_star
    var = sensitivity_results.variances

    #@show mu, var

    results_array = [results_array transpose(mu) transpose(var)]

end

results_array = results_array[:,2:end]

fname = "./src/simulated/sensitivity_results/sensitivity_matrix.dat"
writedlm(fname,results_array)
