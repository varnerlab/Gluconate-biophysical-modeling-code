# in sequential order: delete old files, estimate parameters, plot all figures except for sensitivity, run sensitivity analysis, plot sensitivity array. If not re-estimating parameters, comment out line #16 and #19.

# delete all current simulation files except for the subfolders. 
using Shell
# Shell.run("rm -f src/simulated/gluconate_dynamics/20mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/10mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/5mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/1mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/0.5mM/*")
# Shell.run("rm -f src/simulated/gluconate_dynamics/0.1mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/0.01mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/0.001mM/*");
# Shell.run("rm -f src/simulated/gluconate_dynamics/0.0001mM/*");
# Shell.run("rm -f src/simulated/dose_response_simulations/results_matrix.dat");
# Shell.run("rm -f src/simulated/sensitivity_results/Sensitivity_matrix.dat"); 
# Shell.run("rm -f src/simulated/poets_ensemble/*");

# Estimate parameters
# include("parameter_estimation.jl")

# After estimating parameters
# include("gluconate_dynamics.jl") #MAKE SURE to confirm the simulation file being used from the simulated/poets_ensemble directory, eg. PC_T10.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used from the simulated/poets_ensemble directory, eg. PC_T10.dat
include("plot_protein_venus_automated.jl")
include("plot_protein_gntr.jl")
include("plot_mrna_venus.jl")
include("plot_mrna_gntr.jl")
include("plot_dose_response.jl")

# Sensitivity
include("sensitivity_run.jl")
include("plot_sensitivity_array.jl")

