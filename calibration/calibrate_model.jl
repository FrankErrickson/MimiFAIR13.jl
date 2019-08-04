using Distributions
using Klara

# Load functions needed for MCMC calibration.
include("run_calibration_model.jl")
include("calibration_functions.jl")

# --------------------------------------------------------------------------------------------------
# User-Specified Settings To Change.
# --------------------------------------------------------------------------------------------------

# First year to run model during calibration (note, FAIR defaults to start in 1765 with RCP data).
start_year = 1765

# Final year for model calibration (defaults to 2017).
end_year = 2017

# Number of final samples to draw from joint posterior (after burn-in).
size_chain = 10_000

# Length of burn-in period (number of initial MCMC samples to discard).
burn_in_length = size_chain * 0.1

# Provide a vector of lengths to create additional thinned chain samples (leave empty for no thinning).
# Each sample will be named "thinned_chain_" plus its final thinned length (e.g. "thinned_chain_100000").
thinning = [1000, 100]

# Select initial starting values and step sizes for uncertain parameters (uppre and lower bounds shown for convenience).
parameter_names = ["TCR", "RWF", "F2x_CO₂", "r0", "rC", "rT", "d1", "d2", "aero_direct_rf_scale", "aero_indirect_rf_scale", "CO2_0", "σ_temp", "σ_co2inst", "σ_co2ice", "ρ_temp", "ρ_co2inst"]
parameter_lower = [0.333, 0.3563, 2.81,  29.47, 0.016,  3.507, 113.0, 2.1, -1.11, -1.11, 275.6, 1e-5,  1e-5,  1e-5, 1e-5,  1e-5]
parameter_upper = [7.126, 0.999,  4.61,  40.53, 0.022,  4.823, 365.0, 6.1,  1.33,  1.33, 280.4, 0.2,   1.0,   10.0, 0.999, 0.999]
# Use FAIR reported best fit for TCR and ECS (implied by RWF). Otherwise set as mean of prior distribution.
parameter_start = [1.53,  0.535,  3.71,  35.0,  0.019,  4.165, 239.0, 4.1,  0.11,  0.11, 278.0, 0.1,   0.5,   5.0,  0.5,   0.5]
# Initial step size, crudely ≈ st.dev / 50 for first pass.
mcmc_step       = [0.1,   0.01,   0.025, 0.15,  0.0001, 0.02,  3.5,   0.05, 0.03,  0.03, 0.07,  0.003, 0.015, 0.15, 0.015, 0.015]
# Combine all the parameter settings.
parameter_info = DataFrame(names=parameter_names, lower_bound=parameter_lower, upper_bound=parameter_upper, initial_condition = parameter_start, step_size = mcmc_step)


# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# Run MCMC Calibraiton.
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Create function to run calibration version of FAIR and log-posterior function.
run_fair! = construct_run_calibration_model(start_year=start_year, end_year=end_year)
fair_log_posterior  = construct_log_posterior(run_fair!, start_year=start_year, end_year=end_year)

# Run the MCMC calibration.
mcmc_chain = run_mcmc(size_chain, Int64(burn_in_length), fair_log_posterior, parameter_info[:initial_condition], parameter_info[:step_size])

# If there are user supplied thinning values, create thinned chain samples.
if !isempty(thinning)
	for i = 1:length(thinning)
		# Calculate equally spaced indices for each thinning length.
		thin_indices = trunc(Int64, collect(linspace(1, size(mcmc_chain,2), thinning[i])))
		# Name each variable based on thinned length and then index into full length chain.
		@eval $(Symbol("thinned_chain_"*string(thinning[i]))) = $(mcmc_chain[:,thin_indices])
	end
end



# --------------------------------------------------------------------------------------------------
# Save Results
# --------------------------------------------------------------------------------------------------

# TODO
