# --------------------------------------------------
# Direct radiative forcing effect from aerosols.
# --------------------------------------------------

@defcomp aerosol_indirect_rf begin

    ϕ                     = Parameter()             # Scale factor.
    b_SOx                 = Parameter()             # Sensitivity to sulfur oxides emissions.
    b_POM                 = Parameter()             # Sensitivity to black carbon + organic carbon emissions.
    SOx_emiss_1765        = Parameter()             # 1765 sulfur oxides emissions (Mt yr⁻¹).
    BC_OC_emiss_1765      = Parameter()             # Black carbon + organic carbon emissions (Mt yr⁻¹).
    F_1765                = Parameter()             # For AR5 scaling, pre-industrial forcing was not zero because there were some emissions (use estimates from Skeie et al).
    F_2011                = Parameter()             # 2011 forcing for AR5 scaling (use estimates from Skeie et al).
    scale_AR5             = Parameter{Bool}()       # Scale the forcing output so that the best estimate forcing in 2011 is -0.45 Wm⁻² based on 2011 emissions from the RCPs.
    fix_pre1850_RCP_ai    = Parameter{Bool}()       # Use different relationship for 1750/65 to 1850 based on anthropogenic emissions from Skeie et al (2011) for 1750 (atmos-chem-phys.net/11/11827/2011).
    SOx_emissions         = Parameter(index=[time]) # Sulfur oxides emissions (MtS yr⁻¹).
    BC_emissions          = Parameter(index=[time]) # Black carbon emissions (Mt yr⁻¹).
    OC_emissions          = Parameter(index=[time]) # Organic carbon emissions (Mt yr⁻¹).

    interp_SOx            = Variable(index=[time])  # Linearly interpolated emissions for 1765-1850 (Mt yr⁻¹).
    interp_BC_OC          = Variable(index=[time])  # Linearly interpolated emissions for 1765-1850 (Mt yr⁻¹).
    ERF_aero_cloud        = Variable(index=[time])  # Indirect radiative forcing from aerosols (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Emulation of the global aerosol model of Ghan et al. (2013) to estimate indirect aerosol forcing from precursor emissions.
        if gettime(t) >= 1850 || p.fix_pre1850_RCP_ai == false
            v.ERF_aero_cloud[t] = p.ϕ * log(1.0 + p.b_SOx * p.SOx_emissions[t] + p.b_POM * (p.BC_emissions[t] + p.OC_emissions[t]))
        else
            # Linearly interpolate between 1765 and 1850 if using different relationship for 1750/65 to 1850.
            v.interp_SOx[t] = (gettime(t) - 1765) / 85.0 * p.SOx_emissions[TimestepValue(1850)] + (1850 - gettime(t)) / 85.0 * p.SOx_emiss_1765
            v.interp_BC_OC[t] = (gettime(t) - 1765) / 85.0 * (p.BC_emissions[TimestepValue(1850)] + p.OC_emissions[TimestepValue(1850)]) + (1850 - gettime(t)) / 85.0 * p.BC_OC_emiss_1765

            # Calculate forcing for 1765-1850 using interpolated emission values.
            v.ERF_aero_cloud[t] = p.ϕ * log(1.0 + p.b_SOx * v.interp_SOx[t] + p.b_POM * v.interp_BC_OC[t])
        end

        # From FAIR: "If True, scale the forcing output so that the best estimate forcing in 2011 is -0.45 W/m2 based on 2011 emissions from the RCPs. The Ghan emulator is built on results 
        #             from the CAM5 GCM. As reported in AR5 WG1 Ch7, GCMs tend to overestimate forcing from aerosol-cloud interactions."
        if p.scale_AR5
            v.ERF_aero_cloud[t] = (v.ERF_aero_cloud[t] - p.F_1765) * (-0.45/(p.F_2011 - p.F_1765))
        else
            v.ERF_aero_cloud[t] = (v.ERF_aero_cloud[t] - p.F_1765)
        end
    end
end
