# --------------------------------------------------
# Land use forcing from surface albedo change.
# --------------------------------------------------

@defcomp landuse_rf begin

    landuse_emiss    = Parameter(index=[time]) # RCP 'OtherCO2' emissions (GtC yr⁻¹).

    cumulative_emiss = Variable(index=[time])  # Cumulative 'OtherCO2' emissions (GtC).
    forcing_landuse  = Variable(index=[time])  # Forcing from land use change (Wm⁻²).
end


function run_timestep(s::landuse_rf, t::Int)
    v = s.Variables
    p = s.Parameters

    # Calculate cumulative land use emissions (based on RCP 'OtherCO2' emissions).
    if t == 1
        v.cumulative_emiss[t] = p.landuse_emiss[t]
    else
        v.cumulative_emiss[t] = v.cumulative_emiss[t-1] + p.landuse_emiss[t]
    end

    # Caluclate land use forcing (based on regression of non-fossil CO₂ emissions against AR5 land use forcing).
    v.forcing_landuse[t] = -0.00113789 * v.cumulative_emiss[t]
end
