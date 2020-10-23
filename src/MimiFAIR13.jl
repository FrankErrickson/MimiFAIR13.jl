# Create a module.
module MimiFAIR13

# Load required packages.
using CSVFiles, DataFrames, Mimi, NLsolve

# Load helper functions and FAIR v1.3 components.
include("helper_functions.jl")
include("components/ch4_cycle.jl")
include("components/n2o_cycle.jl")
include("components/co2_cycle.jl")
include("components/other_ghg_cycles.jl")
include("components/ch4_rf.jl")
include("components/n2o_rf.jl")
include("components/other_ghg_rf.jl")
include("components/co2_rf.jl")
include("components/trop_o3_rf.jl")
include("components/strat_o3_rf.jl")
include("components/aerosol_direct_rf.jl")
include("components/aerosol_indirect_rf.jl")
include("components/bc_snow_rf.jl")
include("components/landuse_rf.jl")
include("components/contrails_rf.jl")
include("components/total_rf.jl")
include("components/temperature.jl")

# Export the following functions.
export load_fair_data, create_fair


# -----------------------------------------------------
# Create function to build Mimi FAIR v1.3
# -----------------------------------------------------

function create_fair(;rcp_scenario::String="RCP85", start_year::Int64=1765, end_year::Int64=2500, F2x::Float64=3.71, TCR::Float64=1.6, ECS::Float64=2.75, d::Array{Float64,1}=[239.0, 4.1])

    # ---------------------------------------------
    # Set up some data.
    # ---------------------------------------------

    # Load RCP and other data needed to construct FAIR.
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = load_fair_data(start_year, end_year, rcp_scenario)

    # Names of minor greenhouse gases and ozone-depleting substances.
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
    ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    # Calculate CO₂ radiative forcing scaling term (to keep forcing from 2x CO₂ consistent).
    scale_co2_forcing = co2_rf_scale(F2x, gas_data[findfirst((in)(["CO2"]), gas_data[!, :gas]), :pi_conc], gas_data[findfirst((in)(["N2O"]), gas_data[!, :gas]), :pi_conc])

    # Calcualte coefficients of slow and fast temperature change in each timestep based on user=specified parameters.
    q = calculate_q(TCR, ECS, d, F2x)

    # ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------

    # Create a Mimi model.
    fair = Model()

    #Set number of timesteps.
    nsteps = length(start_year:end_year)

    # Set time index
    set_dimension!(fair, :time, start_year:end_year)

    # Set index for Kyoto and ozone ozone-depleting gases.
    set_dimension!(fair, :other_ghg, other_ghg_names)
    set_dimension!(fair, :ozone_depleting_substances, ods_names)

    # ---------------------------------------------
    # Add components to model.
    # ---------------------------------------------

    add_comp!(fair, ch4_cycle)
    add_comp!(fair, n2o_cycle)
    add_comp!(fair, other_ghg_cycles)
    add_comp!(fair, co2_cycle)
    add_comp!(fair, ch4_rf)
    add_comp!(fair, n2o_rf)
    add_comp!(fair, other_ghg_rf)
    add_comp!(fair, co2_rf)
    add_comp!(fair, trop_o3_rf)
    add_comp!(fair, strat_o3_rf)
    add_comp!(fair, aerosol_direct_rf)
    add_comp!(fair, aerosol_indirect_rf)
    add_comp!(fair, bc_snow_rf)
    add_comp!(fair, landuse_rf)
    add_comp!(fair, contrails_rf)
    add_comp!(fair, total_rf)
    add_comp!(fair, temperature)

    # ---------------------------------------------
    # Assign values to model parameters.
    # ---------------------------------------------

    # ---- Parameters Common to Multiple Components ---- #
    set_param!(fair, :mol_weight_N,    gas_data[findfirst(gas_data[!,:gas] .== "N"), :mol_weight])
    set_param!(fair, :CO₂_0,           gas_data[findfirst(gas_data[!,:gas] .== "CO2"), :pi_conc])
    set_param!(fair, :N₂O_0,           gas_data[findfirst(gas_data[!,:gas] .== "N2O"), :pi_conc])
    set_param!(fair, :CH₄_0,           gas_data[findfirst(gas_data[!,:gas] .== "CH4"), :pi_conc])
    set_param!(fair, :other_ghg_0,     gas_data[findall((in)(other_ghg_names), gas_data[!,:gas]), :pi_conc])
    set_param!(fair, :SOx_emissions,   rcp_emissions.SOx)
    set_param!(fair, :CO_emissions,    rcp_emissions.CO)
    set_param!(fair, :NMVOC_emissions, rcp_emissions.NMVOC)
    set_param!(fair, :BC_emissions,    rcp_emissions.BC)
    set_param!(fair, :NOx_emissions,   rcp_emissions.NOx)
    set_param!(fair, :OC_emissions,    rcp_emissions.OC)

    # ---- Methane Cycle ---- #
    set_param!(fair, :ch4_cycle, :CH₄_fossil_emissions, rcp_emissions.CH4)
    set_param!(fair, :ch4_cycle, :CH₄_natural_emissions, rcp_emissions.NaturalCH4)
    set_param!(fair, :ch4_cycle, :CH₄_τ, 9.3)
    set_param!(fair, :ch4_cycle, :fossil_frac, gas_fractions.ch4_fossil)
    set_param!(fair, :ch4_cycle, :oxidation_frac, 0.61)
    set_param!(fair, :ch4_cycle, :mol_weight_CH₄, gas_data[findfirst(gas_data[!,:gas] .== "CH4"), :mol_weight])
    set_param!(fair, :ch4_cycle, :mol_weight_CO₂, gas_data[findfirst(gas_data[!,:gas] .== "CO2"), :mol_weight])
    set_param!(fair, :ch4_cycle, :emiss2conc_ch4, conversions[findfirst(conversions[!,:gases] .== "CH4"), :emiss2conc])

    # ---- Nitrous Oxide Cycle ---- #
    set_param!(fair, :n2o_cycle, :N₂O_fossil_emissions, rcp_emissions.N2O)
    set_param!(fair, :n2o_cycle, :N₂O_natural_emissions, rcp_emissions.NaturalN2O)
    set_param!(fair, :n2o_cycle, :N₂O_τ, 121.0)
    set_param!(fair, :n2o_cycle, :emiss2conc_n2o, conversions[findfirst(conversions[!,:gases] .== "N2O"), :emiss2conc])


    # ---- Carbon Cycle ---- #
    set_param!(fair, :co2_cycle, :Cacc_0, 0.0)
    set_param!(fair, :co2_cycle, :r0, 35.0)
    set_param!(fair, :co2_cycle, :rC, 0.019)
    set_param!(fair, :co2_cycle, :rT, 4.165)
    set_param!(fair, :co2_cycle, :iIRF_max, 97.0)
    set_param!(fair, :co2_cycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    set_param!(fair, :co2_cycle, :CO₂_τ, [10.0^6, 394.4, 36.54, 4.304])
    set_param!(fair, :co2_cycle, :E, (rcp_emissions[!,:FossilCO2] .+ rcp_emissions[!,:OtherCO2]))
    set_param!(fair, :co2_cycle, :gtc2ppm, conversions[findfirst(conversions[!,:gases] .== "CO2"), :emiss2conc])

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    set_param!(fair, :other_ghg_cycles, :other_ghg_τ, gas_data[findall((in)(other_ghg_names), gas_data[!,:gas]), :lifetimes])
    set_param!(fair, :other_ghg_cycles, :emissions_other_ghg, Matrix(rcp_emissions[:, Symbol.(other_ghg_names)]))
    set_param!(fair, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions[!,:gases]), :emiss2conc])

    # ---- Methane Radiative Forcing ---- #
    set_param!(fair, :ch4_rf, :H₂O_share, 0.12)
    set_param!(fair, :ch4_rf, :scale_CH₄, 1.0)
    set_param!(fair, :ch4_rf, :a₃, -1.3e-6)
    set_param!(fair, :ch4_rf, :b₃, -8.2e-6)

    # ---- Nitrous Oxide Radiative Forcing ---- #
    set_param!(fair, :n2o_rf, :a₂, -8.0e-6)
    set_param!(fair, :n2o_rf, :b₂, 4.2e-6)
    set_param!(fair, :n2o_rf, :c₂, -4.9e-6)

    # ---- Carbon Dioxide Radiative Forcing ---- #
    set_param!(fair, :co2_rf, :a₁, -2.4e-7)
    set_param!(fair, :co2_rf, :b₁, 7.2e-4)
    set_param!(fair, :co2_rf, :c₁, -2.1e-4)
    set_param!(fair, :co2_rf, :scale_CO₂, scale_co2_forcing)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    set_param!(fair, :other_ghg_rf, :radiative_efficiency, gas_data[findall((in)(other_ghg_names), gas_data[!,:gas]), :rad_eff])

    # ---- Tropospheric Ozone Radiative Forcing ---- #
    set_param!(fair, :trop_o3_rf, :mol_weight_NO, gas_data[findfirst(gas_data[!,:gas] .== "NO"), :mol_weight])
    set_param!(fair, :trop_o3_rf, :T0, 0.0)
    set_param!(fair, :trop_o3_rf, :fix_pre1850_RCP_o3, true)

    # ---- Stratospheric Ozone Radiative Forcing ---- #
    set_param!(fair, :strat_o3_rf, :Br, gas_data[findall((in)(ods_names), gas_data[!,:gas]), :br_atoms])
    set_param!(fair, :strat_o3_rf, :Cl, gas_data[findall((in)(ods_names), gas_data[!,:gas]), :cl_atoms])
    set_param!(fair, :strat_o3_rf, :FC, gas_data[findall((in)(ods_names), gas_data[!,:gas]), :strat_frac])
    set_param!(fair, :strat_o3_rf, :δ1, -1.46030698e-5)
    set_param!(fair, :strat_o3_rf, :δ2, 2.05401270e-3)
    set_param!(fair, :strat_o3_rf, :δ3, 1.03143308)
    set_param!(fair, :strat_o3_rf, :ODS₀, gas_data[findall((in)(ods_names), gas_data[!,:gas]), :pi_conc])

    # ---- Aerosol Direct Radiative Forcing ---- #
    set_param!(fair, :aerosol_direct_rf, :β_SOx, -6.2227e-3)
    set_param!(fair, :aerosol_direct_rf, :β_CO, 0.0)
    set_param!(fair, :aerosol_direct_rf, :β_NMVOC, -3.8392e-4)
    set_param!(fair, :aerosol_direct_rf, :β_NOx, -1.16551e-3)
    set_param!(fair, :aerosol_direct_rf, :β_BC, 1.601537e-2)
    set_param!(fair, :aerosol_direct_rf, :β_OC, -1.45339e-3)
    set_param!(fair, :aerosol_direct_rf, :β_NH3, -1.55605e-3)
    set_param!(fair, :aerosol_direct_rf, :NH3_emissions, rcp_emissions.NH3)

    # ---- Aerosol Indirect Radiative Forcing ---- #
    set_param!(fair, :aerosol_indirect_rf, :ϕ, -1.95011431)
    set_param!(fair, :aerosol_indirect_rf, :b_SOx, 0.01107147)
    set_param!(fair, :aerosol_indirect_rf, :b_POM, 0.01387492)
    set_param!(fair, :aerosol_indirect_rf, :SOx_emiss_1765, 1.0)
    set_param!(fair, :aerosol_indirect_rf, :BC_OC_emiss_1765, 11.2)
    set_param!(fair, :aerosol_indirect_rf, :scale_AR5, true)
    set_param!(fair, :aerosol_indirect_rf, :fix_pre1850_RCP_ai, true)
    set_param!(fair, :aerosol_indirect_rf, :F_1765, -0.3002836449793625)
    set_param!(fair, :aerosol_indirect_rf, :F_2011, -1.5236182344467388)

    # ---- Land Use Change Radiative Forcing ---- #
    set_param!(fair, :landuse_rf, :landuse_emissions, rcp_emissions.OtherCO2)

    # ---- Contrails Radiative Forcing ---- #
    set_param!(fair, :contrails_rf, :frac, gas_fractions.nox_aviation)
    set_param!(fair, :contrails_rf, :E_ref, 2.946)
    set_param!(fair, :contrails_rf, :F_ref, 0.0448)
    set_param!(fair, :contrails_rf, :ref_is_NO2, true)
    set_param!(fair, :contrails_rf, :mol_weight_NO₂, gas_data[findfirst(gas_data[!,:gas] .== "NO2"), :mol_weight])

    # ---- Total Radiative Forcing ---- #
    set_param!(fair, :total_rf, :F_volcanic, volcano_forcing)
    set_param!(fair, :total_rf, :F_solar, solar_forcing)
    set_param!(fair, :total_rf, :F_exogenous, zeros(nsteps))
    set_param!(fair, :total_rf, :efficacy_CO₂, 1.0)
    set_param!(fair, :total_rf, :efficacy_CH₄, 1.0)
    set_param!(fair, :total_rf, :efficacy_CH₄_H₂O, 1.0)
    set_param!(fair, :total_rf, :efficacy_N₂O, 1.0)
    set_param!(fair, :total_rf, :efficacy_other_ghg, ones(length(other_ghg_names)))
    set_param!(fair, :total_rf, :efficacy_trop_O₃, 1.0)
    set_param!(fair, :total_rf, :efficacy_strat_O₃, 1.0)
    set_param!(fair, :total_rf, :efficacy_aerosol_direct, 1.0)
    set_param!(fair, :total_rf, :efficacy_aerosol_indirect, 1.0)
    set_param!(fair, :total_rf, :efficacy_bcsnow, 3.0)
    set_param!(fair, :total_rf, :efficacy_landuse, 1.0)
    set_param!(fair, :total_rf, :efficacy_contrails, 0.0) # Note: Efficacy set to 0.0 to match default settings in Python version of FAIR.

    # ---- Global Temperature Anomaly ---- #
    set_param!(fair, :temperature, :d, d)
    set_param!(fair, :temperature, :q, q)
    set_param!(fair, :temperature, :F2x, F2x)

    # ---------------------------------------------
    # Create connections between Mimi components.
    # ---------------------------------------------
    connect_param!(fair, :co2_cycle, :T, :temperature, :T)
    connect_param!(fair, :trop_o3_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(fair, :strat_o3_rf, :conc_ODS, :other_ghg_cycles, :conc_ods)
    connect_param!(fair, :ch4_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(fair, :ch4_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(fair, :n2o_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(fair, :n2o_rf, :CH₄, :ch4_cycle, :CH₄)
    connect_param!(fair, :n2o_rf, :CO₂, :co2_cycle, :C)
    connect_param!(fair, :co2_rf, :N₂O, :n2o_cycle, :N₂O)
    connect_param!(fair, :co2_rf, :CO₂, :co2_cycle, :C)
    connect_param!(fair, :other_ghg_rf, :conc_other_ghg, :other_ghg_cycles, :conc_other_ghg)
    connect_param!(fair, :trop_o3_rf, :temperature, :temperature, :T)
    connect_param!(fair, :total_rf, :F_CO₂, :co2_rf, :rf_co2)
    connect_param!(fair, :total_rf, :F_CH₄, :ch4_rf, :forcing_CH₄)
    connect_param!(fair, :total_rf, :F_CH₄_H₂O, :ch4_rf, :forcing_CH₄_H₂O)
    connect_param!(fair, :total_rf, :F_N₂O, :n2o_rf, :forcing_N₂O)
    connect_param!(fair, :total_rf, :F_other_ghg, :other_ghg_rf, :other_ghg_rf)
    connect_param!(fair, :total_rf, :F_trop_O₃, :trop_o3_rf, :forcing_trop_O₃)
    connect_param!(fair, :total_rf, :F_strat_O₃, :strat_o3_rf, :forcing_strat_O₃)
    connect_param!(fair, :total_rf, :F_aerosol_direct, :aerosol_direct_rf, :F_aerosol_direct)
    connect_param!(fair, :total_rf, :F_aerosol_indirect, :aerosol_indirect_rf, :ERF_aero_cloud)
    connect_param!(fair, :total_rf, :F_bcsnow, :bc_snow_rf, :forcing_BC_snow)
    connect_param!(fair, :total_rf, :F_landuse, :landuse_rf, :forcing_landuse)
    connect_param!(fair, :total_rf, :F_contrails, :contrails_rf, :forcing_contrails)
    connect_param!(fair, :temperature, :F, :total_rf, :total_forcing)

    # Return FAIR model.
    return fair
end

end # module
