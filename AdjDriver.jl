# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

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

# include -
include("Include.jl")

# Setup the timescale of the simulation -
time_start = 0.0
time_stop = 24.0
time_step_size = 0.1
number_of_timesteps = length(time_start:time_step_size:time_stop)

# Load the data dictionary (default parameter/biophysical values) -
data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size))

# What is the size of the system?
number_of_states = data_dictionary["number_of_states"]

# main loop -
parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
average_scaled_sensitivity_array = zeros(number_of_states,1)

# setup progress meter -
number_of_parameters = length(parameter_name_mapping_array)
p = Progress(number_of_parameters,color=:yellow)

for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)

  # grab the dictionary -
  local_data_dictionary = deepcopy(data_dictionary)

  # Solve the adj simulation -
  # You need to point this to your specific adj simulation code -
  #
  # e.g.,
  # (T,X) = adj_washout_inducer(time_start,time_stop,time_step_size,parameter_index,local_data_dictionary)

  # dump the raw sensitivity arrays to disk -
  # you can modify this to point to some place on disk ...
  data_array = [T X]
  file_path = "./sensitivity/AdjSimulation-P"*string(parameter_index)*".dat"
  writedlm(file_path,data_array)

  # message -
  msg = "Completed parameter $(parameter_index) of $(number_of_parameters)"
  ProgressMeter.next!(p; showvalues = [(:status,msg)])
end
