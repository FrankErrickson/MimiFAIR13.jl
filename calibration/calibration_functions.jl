# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# This file contains functions and code to calculate the joint posterior parameter distribution.
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------



#######################################################################################################################
# LOAD CALIBRATION DATA
########################################################################################################################
# Description: This function loads and combines all of the calibration data and observation measurement errors.
#
# Function Arguments:
#
#       start_year:      First year to run the model for calibration (default 1765).
#       end_year:        Final year to run the model for calibration (default 2017).
#----------------------------------------------------------------------------------------------------------------------

function load_calibration_data(;start_year::Int64=1765, end_year::Int64=2017)

    # Create column of calibration years and calculate indicies for calibration time period.
    calibration_data = DataFrame(year = collect(start_year:end_year))
    index_start, index_end = findin(collect(1765:2500), [start_year, end_year])

    #-------------------------------------------------------------------
    # HadCRUT4 temperature data (anomalies relative to 1861-1880 mean).
    #-------------------------------------------------------------------
    raw_temp_data = DataFrame(load(joinpath(@__DIR__, "calibration_data", "hadcrut4_global_temperature.csv"), skiplines_begin=24))

    # Find indices to normalize temperature data to 1861-1880 mean.
    index_1861, index_1880 = findin(raw_temp_data[:year], [1861, 1880])

    # Normalize temperature data to 1861-1880 mean.
    norm_temp_data = DataFrame(year=raw_temp_data[:year], hadcrut_temperature_obs = raw_temp_data[:median] .- mean(raw_temp_data[:median][index_1861:index_1880]))

    # Join data on year.
    calibration_data = join(calibration_data, norm_temp_data, on=:year, kind=:outer)

    # Read in HadCRUT4 1σ errors and rename column.
    raw_temp_errors = DataFrame(load(joinpath(@__DIR__, "calibration_data", "hadcrut4_measurement_errors.csv"), skiplines_begin=21))
    rename!(raw_temp_errors, :one_sigma_all => :hadcrut_temperature_sigma)

    # Join data on year
    calibration_data = join(calibration_data, raw_temp_errors[[:year, :hadcrut_temperature_sigma]], on=:year, kind=:outer)

    #--------------------------------------------------------
    # Mauna Loa instrumental atmospheric CO₂ concentrations.
    #--------------------------------------------------------

    # Load Mauna Loa CO₂ observations and errors, and rename columns.
    raw_mauna_loa_co2_data = DataFrame(load(joinpath(@__DIR__, "calibration_data", "atmospheric_co2_mauna_loa.csv"), skiplines_begin=59))
    rename!(raw_mauna_loa_co2_data, :mean => :maunaloa_co2_obs, :unc => :maunaloa_co2_sigma)

    # Join data on year
    calibration_data = join(calibration_data, raw_mauna_loa_co2_data, on=:year, kind=:outer)

    #-----------------------------------------------------
    # Law Dome ice core atmospheric CO₂ concentrations.
    #-----------------------------------------------------

    # Load Law Dome CO₂ observations and errors, and rename columns.
    raw_law_dome_co2_data = DataFrame(load(joinpath(@__DIR__, "calibration_data", "atmospheric_co2_law_dome.csv"), skiplines_begin=4))
    rename!(raw_law_dome_co2_data, :co2_ice => :lawdome_co2_obs, :one_sigma_error => :lawdome_co2_sigma)

    # Join data on year
    calibration_data = join(calibration_data, raw_law_dome_co2_data, on=:year, kind=:outer)

    # Crop data to appropriate calibration years and return.
    return calibration_data[index_start:index_end, :]
end



#######################################################################################################################
# CALCULATE AR(1) LOG-LIKELIHOOD
########################################################################################################################
# Description: This function calculates the AR(1) log-likelihood in terms of the data-model residuls when accounting for
#              time-varying observation errors. It follows "The Effects of Time-Varying Observation Errors on Semi-Empirical
#              Sea-Level Projections" (Ruckert et al., 2017) DOI 10.1007/s10584-016-1858-z.
#
# Function Arguments:
#
#       residuals: A vector of data-model residuals.
#       σ:         AR(1) innovation standard deviation.
#		ρ:	   	   AR(1) autocorrelation term.
#		ϵ:	       A vector of time-varying observation error estimates (from calibration data sets).
#----------------------------------------------------------------------------------------------------------------------

function hetero_logl_ar1(residuals::Array{Float64,1}, σ::Float64, ρ::Float64, ϵ::Array{Union{Float64, Missings.Missing},1})

    # Calculate length of residuals.
    n=length(residuals)

    # Define AR(1) stationary process variance.
    σ_process = σ^2/(1-ρ^2)

    # Initialize AR(1) covariance matrix (just for convenience).
    H = abs.(collect(1:n)' .- collect(1:n))

    # Calculate residual covariance matrix (sum of AR(1) process variance and observation error variances).
    # Note: This follows Supplementary Information Equation (10) in Ruckert et al. (2017).
    cov_matrix = σ_process * ρ .^ H + diagm(ϵ.^2)

    # Return the log-likelihood.
    # NOTE: This sometimes throws a PosDefException(1) (seems to be an issue with Multi-Variate Normal in Distributions package?)
    # TODO: Figure out this error.
    return logpdf(MvNormal(cov_matrix), residuals)
end



#######################################################################################################################
# CALCULATE LOG-PRIOR
########################################################################################################################
# Description: This function calculates the joint log-prior for a given sample of uncertain model and statistical parameters.
#
# Function Arguments:
#
#       p: A vector of uncertain model and statistical process parameters sampled during calibration.
#----------------------------------------------------------------------------------------------------------------------

function total_log_prior(p::Array{Float64,1})

    # Assign each value in sampled parameter vector to an individual name (for convenience).
    TCR                    = p[1]
    RWF                    = p[2]
    F2x_CO₂                = p[3]
    r0                     = p[4]
    rC                     = p[5]
    rT                     = p[6]
    d1                     = p[7]
    d2                     = p[8]
    aero_direct_rf_scale   = p[9]
    aero_indirect_rf_scale = p[10]
    CO2_0                  = p[11]
    σ_temp                 = p[12]
    σ_co2inst              = p[13]
    σ_co2ice               = p[14]
    ρ_temp                 = p[15]
    ρ_co2inst              = p[16]

    #-----------------------------------------
    # Specify prior parameter distributions.
    #-----------------------------------------

    # *NOTE: For the results below (except for TCR and RWF), when papers provide 5-95% CI intervals, it's assumed this applies to a normal
    #        distribution. The range is then re-scaled so it corresponds to ± 2 SDs and then applied to a uniform distribution. This is just
    #        to get a first pass at MCMC with wide enough priors.

    # From Millar et al. (2017): TCR comes from a log-normal distribution with 5-95% CI of 1.0-2.5. Set prior as the min/max values from 100 million samples from this distribution.
    prior_TCR = Uniform(0.333, 7.126)
    # Millar et al. (2017): Assumes a Gaussian distribution bounded between [0,1] with 5-95% CI of 0.45-0.75.
    # For first pass, just set RWF prior so TCR < ECS, and 0.333 ≤ ECS ≤ 20.0 °C. This could yield some non-physical combinations that throw an error (errors will yield a -Inf probability).
    prior_RWF = Uniform(0.3563, 0.999)
    # FAIR v1.3: 5-95% CI of 13% of default paramter following Millar et al. (2017). NOTE: Python and 2017 paper default parameters do not match.
    prior_r0 = Uniform(29.47, 40.53)
    # FAIR v1.3: 5-95% CI of 13% of default paramter following Millar et al. (2017).
    prior_rC = Uniform(0.016, 0.022)
    # FAIR v1.3: 5-95% CI of 13% of default paramter following Millar et al. (2017).
    prior_rT = Uniform(3.507, 4.823)
    # FAIR v1.3: Assumed mean = 239 years, standard deviation = 63 years.
    prior_d1 = Uniform(113.0, 365.0)
    # FAIR v1.3: Assumed mean = 4.1 years, standard deviation = 1.0 years.
    prior_d2 = Uniform(2.1, 6.1)
    # FAIR v1.3: Assumed Guassian, with 5-95% CI of 20% around best estimate of 3.71 Wm⁻².
    prior_F2x_CO₂ = Uniform(2.81, 4.61)
    # FAIR v1.3: 5-95% range for aerosol effective radiative forcing uncertainty is -89% to +111% from IPCC AR5 WG1 - Table 8.6
        # NOTE: Figure out how FAIR uncertainties map to IPCC uncertainties in Table 8.6 (it isn't clear).
    prior_aero_direct_rf_scale = Uniform(-1.11, 1.33)
    # FAIR v1.3 assumes the same uncertainty for direct and indirect ERF uncertainty. Use the same prior, but allow direct and indirect uncertainties to differ.
    prior_aero_indirect_rf_scale = Uniform(-1.11, 1.33)
    # IPCC notes 1750 CO₂ concentration of 278 ± 2 (5-95% interval). IPCC AR5 WG1 - Table 2.1
    prior_CO2_0 = Uniform(275.6, 280.4)
    # Priors for statistical process parameters taken from past SNEASY work.
    prior_σ_temp = Uniform(1e-5, 0.2)
    prior_σ_co2inst = Uniform(1e-5, 1.0)
    prior_σ_co2ice = Uniform(1e-5, 10.0)
    prior_ρ_temp = Uniform(1e-5, 0.999)
    prior_ρ_co2inst = Uniform(1e-5, 0.999)

    #-------------------------
    # Calculate log-prior.
    #-------------------------
    # Assume prior distribution independence, so total prior is sum of individual log-priors.
    log_prior = logpdf(prior_TCR, TCR) + logpdf(prior_RWF, RWF) + logpdf(prior_r0, r0) + logpdf(prior_rC, rC) + logpdf(prior_rT, rT) + logpdf(prior_d1, d1) + logpdf(prior_d2, d2) + logpdf(prior_F2x_CO₂, F2x_CO₂) + logpdf(prior_aero_direct_rf_scale, aero_direct_rf_scale) + logpdf(prior_aero_indirect_rf_scale, aero_indirect_rf_scale) + logpdf(prior_CO2_0, CO2_0) + logpdf(prior_σ_temp, σ_temp) + logpdf(prior_σ_co2inst, σ_co2inst) + logpdf(prior_σ_co2ice, σ_co2ice) + logpdf(prior_ρ_temp, ρ_temp) + logpdf(prior_ρ_temp, ρ_temp)

    return log_prior
end



#######################################################################################################################
# CONSTRUCT LOG-POSTERIOR FUNCTION
########################################################################################################################
# Description: This function creates a function to calculate the joint log-posterior function given user-specified model
#              settings and sample parameter values. It internally creates a function to calculate FAIR's log-likelihood,
#              and then pairs this with the log-prior function.
#
# Function Arguments:
#
#       f_run_model: A function that takes uncertain FAIR parameters as inputs, and returns the model output being calibrated.
#       start_year:  First year to run model during calibration (defaults to 1765 to work with RCP scenarios).
#       end_year:    Final year for model calibration (defaults to 2017).
#----------------------------------------------------------------------------------------------------------------------

function construct_log_posterior(f_run_model; start_year::Int64=1765, end_year::Int64=2017)

    # Create a vector of calibration years and parameter for number of years.
    calibration_years = collect(start_year:end_year)
    n = length(calibration_years)

    # Load calibration data.
    calibration_data = load_calibration_data(start_year=start_year, end_year=end_year)

    # Calculate indices for each year that has an observation in calibration data sets (some data sets are sparse).
    obs_temperature_indices  = find(x-> !ismissing(x), calibration_data[:hadcrut_temperature_obs])
    obs_co2inst_indices      = find(x-> !ismissing(x), calibration_data[:maunaloa_co2_obs])
    obs_co2ice_indices       = find(x-> !ismissing(x), calibration_data[:lawdome_co2_obs])

    # Allocate arrays to calculate data-model residuals.
    temperature_residual  = zeros(length(obs_temperature_indices))
    co2inst_residual      = zeros(length(obs_co2inst_indices))
    # Used to calculate 8 year mean, centered on year of carbon dioxide ice core observation.
    co2ice_mean           = zeros(length(obs_co2ice_indices))

    # Allocate vectors to store model output across entire calibration period.
    modeled_co2 = zeros(n)
    modeled_temperature = zeros(n)

    #---------------------------------------------------------------------------------------------------------------------------------------
    # Create a function to calculate the log-likelihood for the observations, assuming residual independence across calibration data sets.
    # Here, p is a vector of uncertain model and statistical process parameters sampled during calibration.
    #---------------------------------------------------------------------------------------------------------------------------------------

    function total_log_likelihood(p)
        # Assign each value in sampled parameter vector to an individual name (for convenience).
        TCR                    = p[1]
        RWF                    = p[2]
        F2x_CO₂                = p[3]
        r0                     = p[4]
        rC                     = p[5]
        rT                     = p[6]
        d1                     = p[7]
        d2                     = p[8]
        aero_direct_rf_scale   = p[9]
        aero_indirect_rf_scale = p[10]
        CO2_0                  = p[11]
        σ_temp                 = p[12]
        σ_co2inst              = p[13]
        σ_co2ice               = p[14]
        ρ_temp                 = p[15]
        ρ_co2inst              = p[16]

        # Run FAIR given sampled model parameters, and return model output being calibrated.
        f_run_model(modeled_co2, modeled_temperature, TCR, RWF, F2x_CO₂, r0, rC, rT, d1, d2, aero_direct_rf_scale, aero_indirect_rf_scale, CO2_0)

        #-----------------------------------------------------------------------
        # Global Surface Temperature (normalized to 1861-1880 mean) Likelihood.
        #-----------------------------------------------------------------------
        llik_temp = 0.0

        # Calculate temperature residuals.
        for (i, index)=enumerate(obs_temperature_indices)
            temperature_residual[i] = calibration_data[index, :hadcrut_temperature_obs] - modeled_temperature[index]
        end

        # Calculate temperature likelihood.
        llik_temp = hetero_logl_ar1(temperature_residual, σ_temp, ρ_temp, calibration_data[obs_temperature_indices,:hadcrut_temperature_sigma])

        #-----------------------------------------------------------------------
        # Mauna Loa Atmospheric CO₂ Concentration (instrumental) Likelihood.
        #-----------------------------------------------------------------------
        llik_co2inst = 0.0

        # Calculate Mauna Loa CO₂ residuals.
        for (i, index)=enumerate(obs_co2inst_indices)
            co2inst_residual[i] = calibration_data[index, :maunaloa_co2_obs] - modeled_co2[index]
        end

        # Calculate Mauna Loa CO₂ likelihood.
        llik_co2inst = hetero_logl_ar1(co2inst_residual, σ_co2inst, ρ_co2inst, calibration_data[obs_co2inst_indices, :maunaloa_co2_sigma])

        #-----------------------------------------------------------------------
        # Law Dome Atmospheric CO₂ Concentration (ice core) Likelihood.
        #-----------------------------------------------------------------------
        llik_co2ice = 0.0

        # Calculate Law Dome CO₂ residuals.
        for (i, index)=enumerate(obs_co2ice_indices)
            # Assume Law Dome observations map to 8-year modeled CO₂ mean, centered on year of ice core observation.
            co2ice_mean[i] = mean(modeled_co2[index + (-4:3)])
            # Calculate Law Dome CO₂ likelihood (sparse data, so do not use AR(1) likelihood function).
            llik_co2ice = llik_co2ice + logpdf(Normal(co2ice_mean[i], sqrt(σ_co2ice^2 + calibration_data[index, :lawdome_co2_sigma]^2)), calibration_data[index, :lawdome_co2_obs]-co2ice_mean[i])
        end

        #-----------------------------------------------------------------------
        # Total likelihood.
        #-----------------------------------------------------------------------

        # Assume residual independence, so total log-likelihood is the sum of individual log-likelihoods.
        log_likelihood = llik_temp + llik_co2inst + llik_co2ice

        return log_likelihood
    end

    #---------------------------------------------------------------------------------------------------------------------------------------
    # Create a function to calculate the log-posterior for the observations.
    # Here, p is a vector of uncertain model and statistical process parameters sampled during calibration.
    #---------------------------------------------------------------------------------------------------------------------------------------

    # Log-posterior distribution, with p(θ|data) ∝ L(data|θ) × p(θ).
    function log_posterior(p)
        # Use try-catch because this may fail for (i) PosDefException that can occur with multivariate normal in Distributions pacakge 
        # or (ii) non-physical parameter combination samples that do not allow FAIR to solve for the time-varying α term.
        try
            return total_log_prior(p) + total_log_likelihood(p)
        catch
            return -Inf
        end
    end

    # Return the log-posterior function.
    return log_posterior
end



#######################################################################################################################
# RUN MARKOV CHAIN MONTE CARLO (MCMC) CALIBRATION
########################################################################################################################
# Description: This function uses the Robust Adaptive Metropolis (RAM) MCMC algoirthm (Vihola 2012) to calibrate the
#              uncertain model parameters, give a user-supplied function for calculating the posterior distribution.
#
# Function Arguments:
#
#       n_samples:          # Length of the final Markov chain after burn-in period (each "n" is a single sample from the joint posterior parameter distribution).
#       burn_in_size:       # Number of initial samples in the Markov chain to discard.
#       posterior_function: # Function to calculate the MCMC target distribution (here, this is the posterior distribution).
#       starting_point:     # Initial starting point for MCMC algorithm (a vector of values for each uncertain parameter).
#       mcmc_step_sizes:    # Step-size for MCMC algorithm (a vector of step sizes for each uncertain parameter).
#----------------------------------------------------------------------------------------------------------------------

function run_mcmc(n_samples::Int64, burn_in_size::Int64, posterior_function, starting_point::Array{Float64,1}, mcmc_step_sizes::Array{Float64,1})

    # Set up calibration for joint posterior function using RAM algorithm and user-specified settings.
    model = likelihood_model(BasicContMuvParameter(:p, logtarget=posterior_function), false)
    job = BasicMCJob(model, RAM(diagm(mcmc_step_sizes)), BasicMCRange(nsteps=(n_samples+burn_in_size), burnin=burn_in_size), Dict(:p=>starting_point), verbose=false)
    run(job)

    # Return sampled parameter values (each column = sample from joint posterior).
    return output(job).value
end
