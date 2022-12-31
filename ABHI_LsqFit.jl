using LsqFit
using DelimitedFiles
using Plots
# a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
# model(x, p) will accept the full data set as the first argument `x`.
# This means that we need to write our model function so it applies
# the model to the full dataset. We use `@.` to apply the calculations
# across all rows.
# data_full= readdlm("/Users/abhi/Documents/Research/Writing/A Exam/Figures_A_exam/Dose_response.csv",',')
data_full= readdlm("/Users/abhi/Downloads/trial.csv",',')
data_x=data_full[2:end,1]
data_y=data_full[2:end,2]

data_y_norm = (data_y.-minimum(data_y))/(maximum(data_y)-minimum(data_y))
# @. model(x, p) = p[1]*exp(-x*p[2])
@. model(x, p) = (x^p[1])/(p[2]^p[1] + x^p[1])
p0=[0.5,0.5]


# The function applies the per observation function p[1]*exp(-x[i]*p[2]) to the full dataset in x, with i denoting an observation row. We simulate some data and chose our "true" parameters.

# some example data
# xdata: independent variables
# ydata: dependent variable
# xdata = range(0, stop=10, length=20)
xdata = data_x
# ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
# ydata = model(xdata, data_y)
ydata = model(xdata, data_y_norm)
# p0 = [0.5, 0.5]

# fit = curve_fit(model, xdata, ydata, p0)
fit = curve_fit(model, xdata, ydata, p0)
# fit is a composite type (LsqFitResult), with some interesting values:
#	dof(fit): degrees of freedom
	coefficients=coef(fit) # best fit parameters
#	fit.resid: residuals = vector of residuals
#	fit.jacobian: estimated Jacobian at solution

# We can estimate errors on the fit parameters,
# to get standard error of each parameter:
sigma = stderror(fit)
# to get margin of error and confidence interval of each parameter at 5% significance level:
margin_of_error = margin_error(fit, 0.05)
confidence_inter = confidence_interval(fit, 0.05)

#Try plotting the fit
est_n =coefficients[1]
est_kd = coefficients[2]
f(x) = (x^est_n)/(est_kd^est_n + x^est_n)

range_= [0.0001:0.001:20.1...]
len=length(range_)
model_data=zeros(len,1)
for i in 1:len
    model_data[i]=f(range_[i])
end

plot(log10.(range_),model_data);
scatter!(log10.(data_x),data_y_norm)




# # The finite difference method is used above to approximate the Jacobian.
# # Alternatively, a function which calculates it exactly can be supplied instead.
# function jacobian_model(x,p)
#     J = Array{Float64}(undef, length(x), length(p))
#     @. J[:,1] = exp(-x*p[2])     #dmodel/dp[1]
#     @. @views J[:,2] = -x*p[1]*J[:,1] #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
#     J
# end
# fit = curve_fit(model, jacobian_model, xdata, ydata, p0)