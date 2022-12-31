# ----------------------------------------------------------------------------------- #
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2021-04-29T20:40:53.886
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./CellFree.json", host_type::Symbol = :cell_free)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	GntR
		:gene	;	# 2	Venus
		:gene	;	# 3	sigma_70
		:mrna	;	# 4	mRNA_GntR
		:mrna	;	# 5	mRNA_Venus
		:mrna	;	# 6	mRNA_sigma_70
		:protein	;	# 7	protein_GntR
		:protein	;	# 8	protein_Venus
		:protein	;	# 9	protein_sigma_70
	]

	# we need to store the species symbol array for later -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		996.0	;	# 1	GntR
		720.0	;	# 2	Venus
		720.0	;	# 3	sigma_70
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 4	1	mRNA_GntR
		gene_coding_length_array[2]	;	# 5	2	mRNA_Venus
		gene_coding_length_array[3]	;	# 6	3	mRNA_sigma_70
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 7	1	protein_GntR
		round((0.33)*mRNA_coding_length_array[2])	;	# 8	2	protein_Venus
		round((0.33)*mRNA_coding_length_array[3])	;	# 9	3	protein_sigma_70
	]

	# array of intracellular gene copy numbers -
	gene_abundance_array = [
		0.010	;	# (number/cell) 1	GntR
		0.007	;	# (number/cell) 2	Venus
		0.0	;	# (number/cell) 3	sigma_70
	]

	# initial condition array -
	initial_condition_array = [
		gene_abundance_array[1]	;	# 1	GntR
		gene_abundance_array[2]	;	# 2	Venus
		gene_abundance_array[3]	;	# 3	sigma_70
		0.0	;	# 4	mRNA_GntR
		0.0	;	# 5	mRNA_Venus
		0.0	;	# 6	mRNA_sigma_70
		0.0	;	# 7	protein_GntR
		0.0	;	# 8	protein_Venus
		0.035	;	# 9	protein_sigma_70

		# translation capacity -
		100.0	;	# translation capacity
	]

	binding_parameter_dictionary = Dict{String,Float64}()
	binding_parameter_dictionary["n_GntR_sigma_70"] = 1.0
	binding_parameter_dictionary["K_GntR_sigma_70"] = 0.05
	binding_parameter_dictionary["n_Venus_sigma_70"] = 1.0
	binding_parameter_dictionary["K_Venus_sigma_70"] = 0.05
	binding_parameter_dictionary["n_Venus_GntR"] = 1.0
	binding_parameter_dictionary["K_Venus_GntR"] = 0.05

	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_GntR_RNAP"] = 0.001
	control_parameter_dictionary["W_GntR_sigma_70"] = 1.0
	control_parameter_dictionary["W_Venus_RNAP"] = 0.001
	control_parameter_dictionary["W_Venus_sigma_70"] = 1.0
	control_parameter_dictionary["W_Venus_GntR"] = 1.0
	control_parameter_dictionary["W_sigma_70_RNAP"] = 0.001


	# D- Gluconate parameter values
	gluconate_parameter_dictionary = Dict{String,Float64}()
	gluconate_parameter_dictionary["gluconate_concentration"] = 10 # mM fix at 10mM
	gluconate_parameter_dictionary["n_gluconate_GntR"] = 1.0
	gluconate_parameter_dictionary["K_gluconate_GntR"] = 1.0 # mM

	# degradation modifiers -
	degradation_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	Venus
		0.0	;	# 3	sigma_70
		1.0	;	# 4	mRNA_GntR
		1.0	;	# 5	mRNA_Venus
		1.0	;	# 6	mRNA_sigma_70
		1.0	;	# 7	protein_GntR
		1.0	;	# 8	protein_Venus
		1.0	;	# 9	protein_sigma_70
	]


	# time constant modifiers -
	time_constant_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	Venus
		0.0	;	# 3	sigma_70
		1.0	;	# 4	mRNA_GntR
		1.0	;	# 5	mRNA_Venus
		1.0	;	# 6	mRNA_sigma_70
		1.0	;	# 7	protein_GntR
		3.0	;	# 8	protein_Venus
		1.0	;	# 9	protein_sigma_70
	]

	# Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# Parameter name index array -
	parameter_name_mapping_array = [
		"n_GntR_sigma_70"	;	# 1
		"K_GntR_sigma_70"	;	# 2
		"n_Venus_sigma_70"	;	# 3
		"K_Venus_sigma_70"	;	# 4
		"n_Venus_GntR"	;	# 5
		"K_Venus_GntR"	;	# 6
		"W_GntR_RNAP"	;	# 7
		"W_GntR_sigma_70"	;	# 8
		"W_Venus_RNAP"	;	# 9
		"W_Venus_sigma_70"	;	# 10
		"W_Venus_GntR"	;	# 11
		"W_sigma_70_RNAP"	;	# 12
		"rnapII_concentration"	;	# 13
		"ribosome_concentration"	;	# 14
		"degradation_constant_mRNA"	;	# 15
		"degradation_constant_protein"	;	# 16
		"kcat_transcription"	;	# 17
		"kcat_translation"	;	# 18
		"maximum_specific_growth_rate"	;	# 19
		"saturation_constant_transcription"	;	# 20
		"saturation_constant_translation"	;	# 21
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
	data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #

	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 37.0 	# K
	data_dictionary["half_life_translation_capacity"] = 8.0	 # hr
	data_dictionary["base_line_weight"] = (1.0/1.0)

	#ABHI ADDED MORE parameters
	data_dictionary["transcription_capacity_delay"] = 6.0
	data_dictionary["transcription_capacity_slope"] = 0.5

	data_dictionary["translation_capacity_delay"] = 1.0
	data_dictionary["translation_capacity_slope"] = 0.5

	return data_dictionary
end
