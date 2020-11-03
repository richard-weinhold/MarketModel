# -------------------------------------------------------------------------------------------------
# POMATO - Power Market Tool (C) 2020
# Current Version: 0.2.5
# Created by Richard Weinhold and Robert Mieth
# Licensed under LGPL v3
#
# Language: Julia, v1.3, 1.4, 1.5 (throwing deprecation warnings)
# ----------------------------------
#
# This file:
# Definition of the´MarketModel moduel
# Called by market_model.py, reads pre-processed data from /data_temp/julia_files/data/
# Output: Optimization results saved in /data_temp/julia_files/results/*result_folder*
# -------------------------------------------------------------------------------------------------

module MarketModel
using Logging, LoggingExtras
using DataFrames, CSV, JSON, Dates, Base.Threads
using LinearAlgebra, Distributions, SparseArrays
using JuMP

include("data.jl")
include("pomato.jl")
include("result.jl")
include("data_io.jl")
include("data_processing.jl")
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
