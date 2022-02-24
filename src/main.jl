"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Upper level functions to run the models. These are the functions exposed by MarketModel.
"""


function set_logger()
	# global_logger(ConsoleLogger(stdout, Logging.Info))
	if isfile(pwd()*"/logs/MarketModel.log")
		TeeLogger(MinLevelLogger(FileLogger(pwd()*"/logs/MarketModel.log", append=true), Logging.Info),
		          ConsoleLogger(stdout, Logging.Info)) |> global_logger
		@info("Logfile Found, logging to console and logfile.")
	else
		TeeLogger(ConsoleLogger(stdout, Logging.Info)) |> global_logger
		@info("No logfile Found, logging only to console.")
	end
end

function set_global_optimizer(input_optimizer)
	global optimizer = input_optimizer.Optimizer
	global optimizer_package = input_optimizer
end

function split_timeseries_segments(data::Data, segment_length::Int)
	segments = []
	splits = Int(floor(length(data.t)/segment_length))
	if splits == 0
		return [1:length(data.t)]
	end
	for i in 1:splits - 1
		push!(segments, ((i-1)*segment_length + 1):i*segment_length)
	end
	push!(segments, (splits-1)*segment_length + 1:length(data.t))
 	return segments
end


""" 
	run_market_model(data_dir::String, 
					 result_dir::String, 
					 input_optimizer; 
					 return_result::Bool=false,
                 	 redispatch::Bool=false)

Solves an economic dispatch problem for given data in `data_dir::String` using the supplied solver
supplied as `input_optimizer`. Note, input_optimizer has to be the Julia Package itself. Optionally,
the economic dispatch can be redispatched for feasibility in a given network representation using
the optional argument `redispatch::Bool`. 

Results are saved as in .csv files into `result_dir::String`. Per default the function returns
nothing, but the Result struct can be returned with optional argument `return_result::Bool`. 

"""
function run_market_model(data_dir::String, result_dir::String, input_optimizer;
						  return_result::Bool=false, redispatch::Bool=false)
	# set_logger()
	@info("Read Model Data..")

	data = read_model_data(data_dir)
	data.folders["result_dir"] = result_dir*Dates.format(now(), "dmm_HHMM_SSsss")

	create_folder(data.folders["result_dir"])
	if !("split_timeseries" in keys(data.options))
		data.options["split_timeseries"] = false
	end

	if redispatch
		pomato_results = run_market_model_redispatch(data, input_optimizer)
	else
		pomato_results = run_market_model(data, input_optimizer)
	end

	if return_result
		return pomato_results
	end
	@info("Everything Done!")
end

function run_market_model(data::Data, input_optimizer; 
						  save::Bool=true)


	set_global_optimizer(input_optimizer)	
	pomato_results = Dict{String, Result}()
	if data.options["timeseries"]["split"]
		if true
			@info("Set storage regime in simplified model.")
			set_storage_levels!(data)
		end
		data_full = deepcopy(data)
		model_horizon_segments = split_timeseries_segments(data_full, data.options["timeseries"]["market_horizon"])
		for (i, timesteps) in enumerate(model_horizon_segments)
			data = deepcopy(data_full)
			data.t = data.t[timesteps]
			set_model_horizon!(data, i)
			@info("Initializing Market Model for timestep $(data.t[1].name)...")
			pomato_results[data.t[1].name] = market_model(data).result
		end
	else
		pomato_results[data.t[1].name] = market_model(data).result
	end
	if save
		save_result(concat_results(pomato_results), data.folders["result_dir"])
	end
	return pomato_results
end

function run_market_model_redispatch(data::Data, input_optimizer; 
									 save::Bool=true)

	set_global_optimizer(input_optimizer)
	market_result = concat_results(run_market_model(data, input_optimizer, save=false))
	redispatch_results = redispatch_model(market_result, data)
	if save
		for result in keys(redispatch_results)
			save_result(redispatch_results[result], data.folders["result_dir"]*"_"*result)
		end
	end
	@info("Everything Done!")
	return redispatch_results
end
