include("../mimi_fair_v13.jl")

# Function that creates a function to run Mimi-FAIR given sampled parameters, and return the output being calibrated against observations.
function construct_run_calibration_model(;start_year::Int64=1765, end_year::Int64=2017)

	# Create an instance of FAIR for the calibration time period (default = 1765-2017).
    m = construct_fair(start_year=start_year, end_year=end_year)

    # Set indices for normalizing temerature projections to 1861-1880.
    	# Note from FAIR paper: This represents "...a preindustrial state that is relatively free of volcanic eruptions but with a reasonable global coverage of temperature observations."
    index_1861, index_1880 = findin(collect(start_year:end_year), [1861, 1880])

    # Create the function that runs FAIR and returns appropriate outputs during model calibration.
    function run_calibration_model!(
    	modeled_co2::Vector{Float64},
    	modeled_temperature::Vector{Float64},
    	TCR::Float64,
    	RWF::Float64,
    	F2x_CO₂::Float64,
    	r0::Float64,
    	rC::Float64,
    	rT::Float64,
    	d1::Float64,
    	d2::Float64,
    	aero_direct_rf_scale::Float64,
        aero_indirect_rf_scale::Float64,
    	CO2_0::Float64)

    	#---------------------------------------------
    	# Set uncertain parameters
    	#---------------------------------------------
	    setparameter(m, :co2_cycle, :CO2_0, CO2_0)
	    setparameter(m, :co2_cycle, :r0, r0)
    	setparameter(m, :co2_cycle, :rC, rC)
    	setparameter(m, :co2_cycle, :rT, rT)
    	setparameter(m, :n2o_rf, :CO₂_0, CO2_0)
    	setparameter(m, :co2_rf, :CO₂_0, CO2_0)
		# Calculate CO₂ RF scaling to ensure modeled forcing consistent with sampled F2x_CO₂ value. Assume pre-industrial N₂O constant.
    	setparameter(m, :co2_rf, :rf_scale_CO₂, co2_rf_scale(F2x_CO₂, CO2_0, 273.0))
    	setparameter(m, :temperature, :d, [d1,d2])
		# Calcualte coefficients of slow and fast temperature change in each timestep based on user=specified parameters.
			# Note: ECS calculated based on sampled TCR and realized warming fraction (TCR/ECS ratio).
	    setparameter(m, :temperature, :q, calculate_q(TCR, TCR/RWF, [d1,d2], F2x_CO₂))
	    setparameter(m, :temperature, :F2x, F2x_CO₂)
        # Allow for different uncertainty factors for direct and indirect aerosol radiative forcing.
        setparameter(m, :aerosol_direct_rf, :rf_scale_aero_direct, aero_direct_rf_scale)
        setparameter(m, :aerosol_indirect_rf, :rf_scale_aero_indirect, aero_indirect_rf_scale)

        #---------------------------------------------
        # Caluclate model output used for calibration
        #---------------------------------------------
        run(m)

        # Atmospheric carbon dioxide concentrations.
        modeled_co2[:] = m[:co2_cycle, :C]

        # Global surface temperature anomaly, normalized to 1861-1880 mean.
        modeled_temperature[:] = m[:temperature, :T] .- mean(m[:temperature, :T][index_1861:index_1880])

        return
    end

    return run_calibration_model!
end
