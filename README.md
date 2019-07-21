## Julia Implementation of the FAIR v1.3 Climate Model

Replication of the [Finite Amplitude Impulse-Response simple climate-carbon-cycle model (FAIR) version 1.3](https://github.com/OMS-NetZero/FAIR) in Julia using the [Mimi model development framework](https://github.com/mimiframework/Mimi.jl).

### Software Requirements

This code runs on slightly outdated software, specifically [Julia 0.6](https://julialang.org/downloads/oldreleases.html) and Mimi version 0.4.1. This issue will be addressed in the near future.

### Running The Model
The `user_interface.jl` file contains a number of model settings that can be changed. After modifying and saving this file, the following command will run and solve an instance of the Mimi-FAIR model:
`include("user_interface.jl")`. Results from the solved model, `m`, can be accessed as `m[:component_name, :variable/parameter_name]`, for instance `m[:temperature, :T]` will provide a time series of modeled global surface temperature anomalies.

### Model Validation
The model closely matches the original [Python implementation](https://github.com/OMS-NetZero/FAIR) of FAIR. The plots below compare the Julia and Python versions of FAIR v1.3 across the four RCP scenarios. ![alt text](https://github.com/FrankErrickson/mimi_fair_v13/blob/master/Mimi%20vs%20Python%20FAIR.png)
