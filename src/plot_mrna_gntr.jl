include("Include.jl")

function main(path_to_simulation_dir::String, path_to_plot_file::String)
    # what index is prot gntr?
    state_index = 4
    # how many files?
    searchdir(path,key) = filter(x->contains(x,key),readdir(path))
    file_name_array = searchdir(path_to_simulation_dir, ".dat")
    number_of_trials = length(file_name_array)
    # # initialize -> hardcode the dimension for now
    # data_array = zeros(1201,number_of_trials)
    # time_array = zeros(1201,number_of_trials)

	dt = 0.1 # hr
	tEND = convert(Int64,12/dt)
	t_intervals = collect(0:dt:tEND*dt)
	# t_intervals = [0,2,4,6,8,10,12]
	data_array = zeros(length(t_intervals),number_of_trials)

    # read the simulation dir -
    for (file_index,file_name) in enumerate(file_name_array)
        # load -
        file_path = "$(path_to_simulation_dir)/$(file_name)"
        sim_data_array = readdlm(file_path)

        x = sim_data_array[:,state_index+1]*1000
        t = sim_data_array[:,1]

        spline_obj = DataInterpolations.CubicSpline(x,t)
		mRNA_values = spline_obj.(t_intervals)
		data_array[:, file_index] = mRNA_values

        # PyPlot.plot(t,x,color="dimgrey",alpha=0.80,lw=0.5)
    end
    # plot -
	μ = mean(data_array,dims=2)
    σ = std(data_array,dims=2)
    LB = μ .- (2.576/sqrt(1))*σ
    UB = μ .+ (2.576/sqrt(1))*σ

    fill_between(t_intervals, vec(UB), vec(LB), color="powderblue", alpha=0.80)
    # Plot mean -
    PyPlot.plot(t_intervals,μ,"-",color="black",lw=2)

	# load the experimemtal data -
	experimental_data_dictionary = load_experimental_data_dictionary(pwd())
	# plot the experimemtal data -

	TEXP = experimental_data_dictionary["mRNA_data_array"][[1,2,4,5],1]
	DATA = experimental_data_dictionary["mRNA_data_array"][[1,2,4,5],4]
	STD = (1/sqrt(3))*experimental_data_dictionary["mRNA_data_array"][[1,2,4,5],5]
	yerr_array = transpose([STD STD])
	PyPlot.errorbar(TEXP, DATA,yerr=yerr_array,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=8)

    # labels -
    PyPlot.plot(figsize=(5,4))
    PyPlot.xlabel("Time (hr)", fontsize=20)
    PyPlot.ylabel("[GntR mRNA] (nM)", fontsize=20)
    # PyPlot.axis([-0.5,16.5,-50,1600])
    PyPlot.xticks([0,2,4,6,8,10,12], fontsize=22)
    PyPlot.yticks([0:100:600...], fontsize=22)
    PyPlot.tight_layout()
    PyPlot.savefig(path_to_plot_file)
end

concentration = 10
path_to_simulation_dir = "$(pwd())/src/simulated/gluconate_dynamics/$(concentration)mM"
path_to_plot_file = "$(pwd())/src/plots/mRNA-gntr.pdf"
clf()
main(path_to_simulation_dir, path_to_plot_file)
