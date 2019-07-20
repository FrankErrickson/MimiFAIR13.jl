include("helper_functions.jl")

using Mimi
using DataFrames
using NLsolve
using CSVFiles

include("src/ch4_cycle.jl")
include("src/n2o_cycle.jl")
include("src/co2_cycle.jl")
include("src/other_ghg_cycles.jl")
include("src/ch4_rf.jl")
include("src/n2o_rf.jl")
include("src/other_ghg_rf.jl")
include("src/co2_rf.jl")
include("src/trop_o3_rf.jl")
include("src/strat_o3_rf.jl")
include("src/aerosol_direct_rf.jl")
include("src/aerosol_indirect_rf.jl")
include("src/bc_snow_rf.jl")
include("src/landuse_rf.jl")
include("src/contrails_rf.jl")
include("src/total_rf.jl")
include("src/temperature.jl")

function constructfair(;rcp_scenario::String="RCP85", start_year::Int64=1765, end_year::Int64=2500, F2x::Float64=3.71, TCR::Float64=1.6, ECS::Float64=2.75, d::Array{Float64,1}=[239.0, 4.1])

    # ---------------------------------------------
    # Set up some data.
    # ---------------------------------------------

    # Load RCP and other data needed to construct FAIR.
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = load_fair_data(start_year, end_year, rcp_scenario)

    # Names of minor greenhouse gases and ozone-depleting substances.
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
    ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    # Calculate CO₂ radiative forcing scaling term (to keep forcing from 2x CO₂ consistent).
    scale_co2_forcing = co2_rf_scale(F2x, gas_data[findin(gas_data[:gas],["CO2"]), :pi_conc][1], gas_data[findin(gas_data[:gas],["N2O"]), :pi_conc][1])

    # Calcualte coefficients of slow and fast temperature change in each timestep based on user=specified parameters.
    q = calculate_q(TCR, ECS, d, F2x)

    # ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------

    # Create a Mimi model.
    m = Model()

    #Set number of timesteps.
    nsteps = length(start_year:end_year)

    # Set time index
    setindex(m, :time, nsteps)

    # Set index for Kyoto and ozone ozone-depleting gases.
    setindex(m, :other_ghg, other_ghg_names)
    setindex(m, :ozone_depleting_substances, ods_names)

    # ---------------------------------------------
    # Add components to model.
    # ---------------------------------------------

    addcomponent(m, ch4_cycle)
    addcomponent(m, n2o_cycle)
    addcomponent(m, other_ghg_cycles)
    addcomponent(m, co2_cycle)
    addcomponent(m, ch4_rf)
    addcomponent(m, n2o_rf)
    addcomponent(m, other_ghg_rf)
    addcomponent(m, co2_rf)
    addcomponent(m, trop_o3_rf)
    addcomponent(m, strat_o3_rf)
    addcomponent(m, aerosol_direct_rf)
    addcomponent(m, aerosol_indirect_rf)
    addcomponent(m, bc_snow_rf)
    addcomponent(m, landuse_rf)
    addcomponent(m, contrails_rf)
    addcomponent(m, total_rf)
    addcomponent(m, temperature)

    # ---------------------------------------------
    # Assign values to model parameters.
    # ---------------------------------------------

    # ---- Methane Cycle ---- #
    setparameter(m, :ch4_cycle, :fossil_emiss, rcp_emissions[:CH4])
    setparameter(m, :ch4_cycle, :natural_emiss, rcp_emissions[:NaturalCH4])
    setparameter(m, :ch4_cycle, :CH₄_0, gas_data[findin(gas_data[:gas],["CH4"]), :pi_conc][1])
    setparameter(m, :ch4_cycle, :τ, 9.3)
    setparameter(m, :ch4_cycle, :fossil_frac, gas_fractions[:ch4_fossil])
    setparameter(m, :ch4_cycle, :oxidation_frac, 0.61)
    setparameter(m, :ch4_cycle, :mol_weight_CH₄, gas_data[findin(gas_data[:gas],["CH4"]), :mol_weight][1])
    setparameter(m, :ch4_cycle, :mol_weight_CO₂, gas_data[findin(gas_data[:gas],["CO2"]), :mol_weight][1])
    setparameter(m, :ch4_cycle, :emiss2conc_ch4, conversions[findin(conversions[:gases], ["CH4"]), :emiss2conc][1])

    # ---- Nitrous Oxide Cycle ---- #
    setparameter(m, :n2o_cycle, :fossil_emiss, rcp_emissions[:N2O])
    setparameter(m, :n2o_cycle, :natural_emiss, rcp_emissions[:NaturalN2O])
    setparameter(m, :n2o_cycle, :N₂O_0, gas_data[findin(gas_data[:gas],["N2O"]), :pi_conc][1])
    setparameter(m, :n2o_cycle, :τ, 121.0)
    setparameter(m, :n2o_cycle, :emiss2conc_n2o, conversions[findin(conversions[:gases], ["N2O"]), :emiss2conc][1])

    # ---- Carbon Cycle ---- #
    setparameter(m, :co2_cycle, :CO2_0, gas_data[findin(gas_data[:gas],["CO2"]), :pi_conc][1])
    setparameter(m, :co2_cycle, :r0, 35.0)
    setparameter(m, :co2_cycle, :rC, 0.019)
    setparameter(m, :co2_cycle, :rT, 4.165)
    setparameter(m, :co2_cycle, :iIRF_max, 97.0)
    setparameter(m, :co2_cycle, :a, [0.2173, 0.2240, 0.2824, 0.2763])
    setparameter(m, :co2_cycle, :τ, [10.0^6, 394.4, 36.54, 4.304])
    setparameter(m, :co2_cycle, :E, (rcp_emissions[:FossilCO2] .+ rcp_emissions[:OtherCO2]))
    setparameter(m, :co2_cycle, :gtc2ppm, conversions[findin(conversions[:gases], ["CO2"]), :emiss2conc][1])

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    setparameter(m, :other_ghg_cycles, :τ, gas_data[findin(gas_data[:gas], other_ghg_names), :lifetimes])
    setparameter(m, :other_ghg_cycles, :other_ghg_0, gas_data[findin(gas_data[:gas], other_ghg_names), :pi_conc])
    setparameter(m, :other_ghg_cycles, :emiss_other_ghg, Array(rcp_emissions[Symbol.(other_ghg_names)]))
    setparameter(m, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findin(conversions[:gases], other_ghg_names), :emiss2conc])

    # ---- Methane Radiative Forcing ---- #
    setparameter(m, :ch4_rf, :N₂O_0, gas_data[findin(gas_data[:gas],["N2O"]), :pi_conc][1])
    setparameter(m, :ch4_rf, :H₂O_share, 0.12)
    setparameter(m, :ch4_rf, :scale_CH₄, 1.0)
    setparameter(m, :ch4_rf, :CH₄_0, gas_data[findin(gas_data[:gas],["CH4"]), :pi_conc][1])
    setparameter(m, :ch4_rf, :a₃, -1.3e-6)
    setparameter(m, :ch4_rf, :b₃, -8.2e-6)

    # ---- Nitrous Oxide Radiative Forcing ---- #
    setparameter(m, :n2o_rf, :a₂, -8.0e-6)
    setparameter(m, :n2o_rf, :b₂, 4.2e-6)
    setparameter(m, :n2o_rf, :c₂, -4.9e-6)
    setparameter(m, :n2o_rf, :N₂O_0, gas_data[findin(gas_data[:gas],["N2O"]), :pi_conc][1])
    setparameter(m, :n2o_rf, :CO₂_0, gas_data[findin(gas_data[:gas],["CO2"]), :pi_conc][1])
    setparameter(m, :n2o_rf, :CH₄_0, gas_data[findin(gas_data[:gas],["CH4"]), :pi_conc][1])

    # ---- Carbon Dioxide Radiative Forcing ---- #
    setparameter(m, :co2_rf, :a₁, -2.4e-7)
    setparameter(m, :co2_rf, :b₁, 7.2e-4)
    setparameter(m, :co2_rf, :c₁, -2.1e-4)
    setparameter(m, :co2_rf, :CO₂_0, gas_data[findin(gas_data[:gas],["CO2"]), :pi_conc][1])
    setparameter(m, :co2_rf, :N₂O_0, gas_data[findin(gas_data[:gas],["N2O"]), :pi_conc][1])
    setparameter(m, :co2_rf, :scale_CO₂, scale_co2_forcing)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    setparameter(m, :other_ghg_rf, :other_ghg_0,  gas_data[findin(gas_data[:gas], other_ghg_names), :pi_conc])
    setparameter(m, :other_ghg_rf, :radiative_efficiency,  gas_data[findin(gas_data[:gas], other_ghg_names), :rad_eff])

    # ---- Tropospheric Ozone Radiative Forcing ---- #
    setparameter(m, :trop_o3_rf, :NOx_emissions, rcp_emissions[:NOx])
    setparameter(m, :trop_o3_rf, :CO_emissions, rcp_emissions[:CO])
    setparameter(m, :trop_o3_rf, :NMVOC_emissions, rcp_emissions[:NMVOC])
    setparameter(m, :trop_o3_rf, :mol_weight_N, gas_data[findin(gas_data[:gas],["N"]), :mol_weight][1])
    setparameter(m, :trop_o3_rf, :mol_weight_NO, gas_data[findin(gas_data[:gas],["NO"]), :mol_weight][1])
    setparameter(m, :trop_o3_rf, :CH₄_0, gas_data[findin(gas_data[:gas],["CH4"]), :pi_conc][1])
    setparameter(m, :trop_o3_rf, :T0, 0.0)
    setparameter(m, :trop_o3_rf, :model_years, collect(start_year:end_year))
    setparameter(m, :trop_o3_rf, :fix_pre1850_RCP, true)

    # ---- Stratospheric Ozone Radiative Forcing ---- #
    setparameter(m, :strat_o3_rf, :Br, gas_data[findin(gas_data[:gas], ods_names), :br_atoms])
    setparameter(m, :strat_o3_rf, :Cl, gas_data[findin(gas_data[:gas], ods_names), :cl_atoms])
    setparameter(m, :strat_o3_rf, :FC, gas_data[findin(gas_data[:gas], ods_names), :strat_frac])
    setparameter(m, :strat_o3_rf, :δ1, -1.46030698e-5)
    setparameter(m, :strat_o3_rf, :δ2, 2.05401270e-3)
    setparameter(m, :strat_o3_rf, :δ3, 1.03143308)
    setparameter(m, :strat_o3_rf, :ODS₀, gas_data[findin(gas_data[:gas], ods_names), :pi_conc])

    # ---- Aerosol Direct Radiative Forcing ---- #
    setparameter(m, :aerosol_direct_rf, :β_SOx, -6.2227e-3)
    setparameter(m, :aerosol_direct_rf, :β_CO, 0.0)
    setparameter(m, :aerosol_direct_rf, :β_NMVOC, -3.8392e-4)
    setparameter(m, :aerosol_direct_rf, :β_NOx, -1.16551e-3)
    setparameter(m, :aerosol_direct_rf, :β_BC, 1.601537e-2)
    setparameter(m, :aerosol_direct_rf, :β_OC, -1.45339e-3)
    setparameter(m, :aerosol_direct_rf, :β_NH3, -1.55605e-3)
    setparameter(m, :aerosol_direct_rf, :SOx_emiss, rcp_emissions[:SOx])
    setparameter(m, :aerosol_direct_rf, :CO_emiss, rcp_emissions[:CO])
    setparameter(m, :aerosol_direct_rf, :NMVOC_emiss, rcp_emissions[:NMVOC])
    setparameter(m, :aerosol_direct_rf, :NOx_emiss, rcp_emissions[:NOx])
    setparameter(m, :aerosol_direct_rf, :BC_emiss, rcp_emissions[:BC])
    setparameter(m, :aerosol_direct_rf, :OC_emiss, rcp_emissions[:OC])
    setparameter(m, :aerosol_direct_rf, :NH3_emiss, rcp_emissions[:NH3])

    # ---- Aerosol Indirect Radiative Forcing ---- #
    setparameter(m, :aerosol_indirect_rf, :ϕ, -1.95011431)
    setparameter(m, :aerosol_indirect_rf, :b_SOx, 0.01107147)
    setparameter(m, :aerosol_indirect_rf, :b_POM, 0.01387492)
    setparameter(m, :aerosol_indirect_rf, :SOx_emiss, rcp_emissions[:SOx])
    setparameter(m, :aerosol_indirect_rf, :BC_emiss, rcp_emissions[:BC])
    setparameter(m, :aerosol_indirect_rf, :OC_emiss, rcp_emissions[:OC])
    setparameter(m, :aerosol_indirect_rf, :rcp_1850_index, findin(collect((1765:2500)), 1850)[1])
    setparameter(m, :aerosol_indirect_rf, :SOx_emiss_1765, 1.0)
    setparameter(m, :aerosol_indirect_rf, :BC_OC_emiss_1765, 11.2)
    setparameter(m, :aerosol_indirect_rf, :model_years, collect(start_year:end_year))
    setparameter(m, :aerosol_indirect_rf, :scale_AR5, true)
    setparameter(m, :aerosol_indirect_rf, :fix_pre1850_RCP, true)
    setparameter(m, :aerosol_indirect_rf, :F_1765, -0.3002836449793625)
    setparameter(m, :aerosol_indirect_rf, :F_2011, -1.5236182344467388)

    # ---- Black Carbon on Snow Radiative Forcing ---- #
    setparameter(m, :bc_snow_rf, :BC_emiss, rcp_emissions[:BC])

    # ---- Land Use Change Radiative Forcing ---- #
    setparameter(m, :landuse_rf, :landuse_emiss, rcp_emissions[:OtherCO2])

    # ---- Contrails Radiative Forcing ---- #
    setparameter(m, :contrails_rf, :NOx_emiss, rcp_emissions[:NOx])
    setparameter(m, :contrails_rf, :frac, gas_fractions[:nox_aviation])
    setparameter(m, :contrails_rf, :E_ref, 2.946)
    setparameter(m, :contrails_rf, :F_ref, 0.0448)
    setparameter(m, :contrails_rf, :ref_is_NO2, true)
    setparameter(m, :contrails_rf, :mol_weight_NO₂, gas_data[findin(gas_data[:gas],["NO2"]), :mol_weight][1])
    setparameter(m, :contrails_rf, :mol_weight_N, gas_data[findin(gas_data[:gas],["N"]), :mol_weight][1])

    # ---- Total Radiative Forcing ---- #
    setparameter(m, :total_rf, :F_volcanic, volcano_forcing)
    setparameter(m, :total_rf, :F_solar, solar_forcing)
    setparameter(m, :total_rf, :F_exogenous, zeros(nsteps))
    setparameter(m, :total_rf, :efficacy_CO₂, 1.0)
    setparameter(m, :total_rf, :efficacy_CO₂, 1.0)
    setparameter(m, :total_rf, :efficacy_CH₄, 1.0)
    setparameter(m, :total_rf, :efficacy_CH₄_H₂O, 1.0)
    setparameter(m, :total_rf, :efficacy_N₂O, 1.0)
    setparameter(m, :total_rf, :efficacy_other_ghg, ones(length(other_ghg_names)))
    setparameter(m, :total_rf, :efficacy_trop_O₃, 1.0)
    setparameter(m, :total_rf, :efficacy_strat_O₃, 1.0)
    setparameter(m, :total_rf, :efficacy_aerosol_direct, 1.0)
    setparameter(m, :total_rf, :efficacy_aerosol_indirect, 1.0)
    setparameter(m, :total_rf, :efficacy_bcsnow, 3.0)
    setparameter(m, :total_rf, :efficacy_landuse, 1.0)
    setparameter(m, :total_rf, :efficacy_contrails, 0.0) # Note: Efficacy set to 0.0 to match default settings in Python version of FAIR.

    # ---- Global Temperature Anomaly ---- #
    setparameter(m, :temperature, :d, d)
    setparameter(m, :temperature, :q, q)
    setparameter(m, :temperature, :F2x, F2x)

    # ---------------------------------------------
    # Create connections between Mimi components.
    # ---------------------------------------------
    connectparameter(m, :co2_cycle, :T, :temperature, :T)
    connectparameter(m, :trop_o3_rf, :CH₄, :ch4_cycle, :CH₄)
    connectparameter(m, :strat_o3_rf, :conc_ODS, :other_ghg_cycles, :conc_ods)
    connectparameter(m, :ch4_rf, :CH₄, :ch4_cycle, :CH₄)
    connectparameter(m, :ch4_rf, :N₂O, :n2o_cycle, :N₂O)
    connectparameter(m, :n2o_rf, :N₂O, :n2o_cycle, :N₂O)
    connectparameter(m, :n2o_rf, :CH₄, :ch4_cycle, :CH₄)
    connectparameter(m, :n2o_rf, :CO₂, :co2_cycle, :C)
    connectparameter(m, :co2_rf, :N₂O, :n2o_cycle, :N₂O)
    connectparameter(m, :co2_rf, :CO₂, :co2_cycle, :C)
    connectparameter(m, :other_ghg_rf, :conc_other_ghg, :other_ghg_cycles, :conc_other_ghg)
    connectparameter(m, :trop_o3_rf, :temperature, :temperature, :T)
    connectparameter(m, :total_rf, :F_CO₂, :co2_rf, :rf_co2)
    connectparameter(m, :total_rf, :F_CH₄, :ch4_rf, :forcing_CH₄)
    connectparameter(m, :total_rf, :F_CH₄_H₂O, :ch4_rf, :forcing_CH₄_H₂O)
    connectparameter(m, :total_rf, :F_N₂O, :n2o_rf, :forcing_N₂O)
    connectparameter(m, :total_rf, :F_other_ghg, :other_ghg_rf, :other_ghg_rf)
    connectparameter(m, :total_rf, :F_trop_O₃, :trop_o3_rf, :forcing_trop_O₃)
    connectparameter(m, :total_rf, :F_strat_O₃, :strat_o3_rf, :forcing_strat_O₃)
    connectparameter(m, :total_rf, :F_aerosol_direct, :aerosol_direct_rf, :F_aerosol_direct)
    connectparameter(m, :total_rf, :F_aerosol_indirect, :aerosol_indirect_rf, :ERF_aero_cloud)
    connectparameter(m, :total_rf, :F_bcsnow, :bc_snow_rf, :forcing_BC_snow)
    connectparameter(m, :total_rf, :F_landuse, :landuse_rf, :forcing_landuse)
    connectparameter(m, :total_rf, :F_contrails, :contrails_rf, :forcing_contrails)
    connectparameter(m, :temperature, :F, :total_rf, :total_forcing)

    return m
end
