# --------------------------------------------------
# Land use forcing from surface albedo change.
# --------------------------------------------------

@defcomp landuse_rf begin

    landuse_emissions    = Parameter(index=[time]) # RCP 'OtherCO2' emissions (GtC yr⁻¹).

    cumulative_emissions = Variable(index=[time])  # Cumulative 'OtherCO2' emissions (GtC).
    forcing_landuse      = Variable(index=[time])  # Forcing from land use change (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate cumulative land use emissions (based on RCP 'OtherCO2' emissions).
        if is_first(t)
            v.cumulative_emissions[t] = p.landuse_emissions[t]
        else
            v.cumulative_emissions[t] = v.cumulative_emissions[t-1] + p.landuse_emissions[t]
        end

        # Caluclate land use forcing (based on regression of non-fossil CO₂ emissions against AR5 land use forcing).
        v.forcing_landuse[t] = -0.00113789 * v.cumulative_emissions[t]
    end
end
