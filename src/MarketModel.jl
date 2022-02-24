"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Definition of the MarketModel module. 
Inclusion of all components and modules. 
Expose core functionality. 
"""

module MarketModel
using Logging, LoggingExtras
using DataFrames, CSV, JSON, Dates, Base.Threads
using LinearAlgebra, Distributions, SparseArrays
using Clustering, Interpolations, RollingFunctions
using JuMP

include("data.jl")
include("pomato.jl")
include("result.jl")
include("data_io.jl")
include("data_processing.jl")
include("gurobi_functions.jl")
include("storage_model.jl")
include("pomato_models.jl")
include("model_functions.jl")
include("main.jl")

export run_market_model, run_market_model_redispatch

function __init__()
	# global_logger(ConsoleLogger(stdout, Logging.Info))
	# @info("No arguments passed or not running in repl, initializing in pwd()")
	# @info("Initialized")
	println("Initialized")
end
end # module
