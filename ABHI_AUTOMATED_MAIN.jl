# in sequential order: delete old files, estimate parameters, sensitivity. If not estimating parameters, comment out line #14 and #17.

# delete all current simulation and other files
using Shell
Shell.run("rm -f gluconate_dynamics/20mM/*");
Shell.run("rm -f gluconate_dynamics/10mM/*");
Shell.run("rm -f gluconate_dynamics/5mM/*");
Shell.run("rm -f gluconate_dynamics/1mM/*");
Shell.run("rm -f gluconate_dynamics/0.5mM/*")
Shell.run("rm -f gluconate_dynamics/0.1mM/*");
Shell.run("rm -f gluconate_dynamics/0.01mM/*");
Shell.run("rm -f gluconate_dynamics/0.001mM/*");
Shell.run("rm -f gluconate_dynamics/0.0001mM/*");
# Shell.run("rm -f poets_ensemble_W/*");

# Estimate parameters
# include("Parameter_estimation_W_splined.jl")

# After estimating parameters
include("Gluconate_dynamics_automated.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("prot-venus-plot_automated.jl")
include("prot-GNTR-plot.jl")
include("mrna-venus-plot.jl")
include("mrna-gntr-plot.jl")
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("dose_response_plot_edited.jl")

# Sensitivity
# include("sensitivity-gntr-gluconate-updated.jl")
# include("visualize-sensitivity-array.jl")

