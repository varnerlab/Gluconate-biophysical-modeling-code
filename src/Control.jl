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
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-04-29T20:40:54.294
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = zeros(3)

	# Alias the species -
	GntR = x[1]
	Venus = x[2]
	sigma_70 = x[3]
	mRNA_GntR = x[4]
	mRNA_Venus = x[5]
	mRNA_sigma_70 = x[6]
	protein_GntR = x[7]
	protein_Venus = x[8]
	protein_sigma_70 = x[9]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_GntR_sigma_70 = binding_parameter_dictionary["n_GntR_sigma_70"]
	K_GntR_sigma_70 = binding_parameter_dictionary["K_GntR_sigma_70"]
	n_Venus_sigma_70 = binding_parameter_dictionary["n_Venus_sigma_70"]
	K_Venus_sigma_70 = binding_parameter_dictionary["K_Venus_sigma_70"]
	n_Venus_GntR = binding_parameter_dictionary["n_Venus_GntR"]
	K_Venus_GntR = binding_parameter_dictionary["K_Venus_GntR"]

	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_GntR_RNAP = control_parameter_dictionary["W_GntR_RNAP"]
	W_GntR_sigma_70 = control_parameter_dictionary["W_GntR_sigma_70"]
	W_Venus_RNAP = control_parameter_dictionary["W_Venus_RNAP"]
	W_Venus_sigma_70 = control_parameter_dictionary["W_Venus_sigma_70"]
	W_Venus_GntR = control_parameter_dictionary["W_Venus_GntR"]
	W_sigma_70_RNAP = control_parameter_dictionary["W_sigma_70_RNAP"]

	# Gluconate binds to GntR like allolactose to LacI
	# f_bound represents the fraction of bound Gluconate [Gluconate]^n/([Gluconate]^n+K^n)
	# unbound GntR is represented by (1-f_bound)[GntR]
	# Gluconate is passed as an parameter
	gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
	protein_gluconate = gluconate_parameter_dictionary["gluconate_concentration"]
	n_gluconate_GntR = gluconate_parameter_dictionary["n_gluconate_GntR"]
	K_gluconate_GntR = gluconate_parameter_dictionary["K_gluconate_GntR"]

	f_bound = (protein_gluconate^(n_gluconate_GntR))/(protein_gluconate^(n_gluconate_GntR)+K_gluconate_GntR^(n_gluconate_GntR))
	protein_GntR = (1-f_bound)*protein_GntR

	# Transfer function target:GntR actor:sigma_70
	actor_set_GntR_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_GntR_sigma_70)
	# actor = abs(prod(actor_set_GntR_sigma_70))
	b_GntR_sigma_70 = (actor^(n_GntR_sigma_70))/(K_GntR_sigma_70^(n_GntR_sigma_70)+actor^(n_GntR_sigma_70))

	# Control function for GntR -
	control_array[1] = (W_GntR_RNAP+W_GntR_sigma_70*b_GntR_sigma_70)/(1+W_GntR_RNAP+W_GntR_sigma_70*b_GntR_sigma_70)

	# Transfer function target:Venus actor:sigma_70
	actor_set_Venus_sigma_70 = [
		protein_sigma_70
	]
	actor = prod(actor_set_Venus_sigma_70)
	# actor = abs(prod(actor_set_Venus_sigma_70))
	b_Venus_sigma_70 = (actor^(n_Venus_sigma_70))/(K_Venus_sigma_70^(n_Venus_sigma_70)+actor^(n_Venus_sigma_70))

	# Transfer function target:Venus actor:GntR
	actor_set_Venus_GntR = [
		protein_GntR
	]
	actor = prod(actor_set_Venus_GntR)
	# actor = abs(prod(actor_set_Venus_GntR))
	b_Venus_GntR = (actor^(n_Venus_GntR))/(K_Venus_GntR^(n_Venus_GntR)+actor^(n_Venus_GntR))

	# Control function for Venus -
	control_array[2] = (W_Venus_RNAP+W_Venus_sigma_70*b_Venus_sigma_70)/(1+W_Venus_RNAP+W_Venus_sigma_70*b_Venus_sigma_70+W_Venus_GntR*b_Venus_GntR)

	# Control function for sigma_70 -
	control_array[3] = (W_sigma_70_RNAP)/(1+W_sigma_70_RNAP)

	# #manually add txtl resources decay term for now (logistic decay)
	# f(y)=(1+exp((-8*pi)/(0.4*sqrt(3))))/(1+exp(((y-8)*pi)/(0.4*sqrt(3))))
	# correction_term_manual=f(t)
	# control_array=control_array*correction_term_manual

	# estimate transcription_capacity_delay and transcription_capacity_slope parameters
	transcription_capacity_delay = data_dictionary["transcription_capacity_delay"] 
	transcription_capacity_slope = data_dictionary["transcription_capacity_slope"] 
	# build control array with these two logistic function parameters
	f(y)= (1+exp((-transcription_capacity_delay*pi)/(transcription_capacity_slope*sqrt(3))))/(1+exp(((y-transcription_capacity_delay)*pi)/(transcription_capacity_slope*sqrt(3))))
	correction_term_manual=f(t)
	control_array=control_array*correction_term_manual

	
	# correction_term = (x[10]/100.0)
    # control_array = control_array*correction_term

	# return -
	return control_array
end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-04-29T20:40:54.447
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(3)

	# #manually add txtl resources decay term for now (logistic decay)
	# f(y)=(1+exp((-6*pi)/(0.6*sqrt(3))))/(1+exp(((y-6)*pi)/(0.6*sqrt(3))))
	# correction_term_manual=f(t)
	# control_array=control_array*correction_term_manual

	# estimate transcription_capacity_delay and transcription_capacity_slope parameters
	translation_capacity_delay = data_dictionary["translation_capacity_delay"] 
	translation_capacity_slope = data_dictionary["translation_capacity_slope"] 
	# build control array with these two logistic function parameters
	f(y)= (1+exp((-translation_capacity_delay*pi)/(translation_capacity_slope*sqrt(3))))/(1+exp(((y-translation_capacity_delay)*pi)/(translation_capacity_slope*sqrt(3))))
	correction_term_manual=f(t)
	# control_array=control_array*correction_term_manual
	
	# # correct for "resource"?
	correction_term = (x[10]/100.0)
    control_array = control_array*correction_term
	# # # return -
	return control_array
end
