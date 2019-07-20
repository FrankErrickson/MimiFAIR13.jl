# Load Mimi FAIR v1.3 source code.
include("mimi_fair_v13.jl")

####################################################################################################
# MIMI-FAIR PARAMETERS TO CHANGE
####################################################################################################

# The first year to run the model.
start_year = 1765

# The last year to run the model.
end_year = 2500

# A string indicating which RCP emissions/forcing scenario to use (options = "RCP26", "RCP45", "RCP60", & "RCP85").
rcp = "RCP85"

# Effective radiative forcing from a doubling of carbon dioxide concentrations (Wm⁻²).
f2x = 3.71

# Transient climate response (K).
tcr = 1.6

# Equilibrium climate sensitivity (K).
ecs = 2.75

# Two-element array of coefficients governing slow and fast thermal response times for upper and deep oceans (years).
d = [239.0, 4.1]



####################################################################################################
# RUN MODEL GIVEN USER PARAMETER SETTINGS
####################################################################################################

# Create instance of Mimi-FAIR model given user settings.
m = constructfair(rcp_scenario=rcp, start_year=start_year, end_year=end_year, F2x=f2x, TCR=tcr, ECS=ecs, d=d)

# Run model.
run(m)
