""" Run the models, defined in create_model.jl"""

function set_logger()
	# global_logger(ConsoleLogger(stdout, Logging.Info))
	if isfile(pwd()*"/logs/RedundancyRemoval.log")
		TeeLogger(MinLevelLogger(FileLogger(pwd()*"/logs/MarketModel.log", append=true), Logging.Info),
		          ConsoleLogger(stdout, Logging.Info)) |> global_logger
		@info("Logfile Found, logging to console and logfile.")
	else
		TeeLogger(ConsoleLogger(stdout, Logging.Info)) |> global_logger
		@info("No logfile Found, logging only to console.")
	end
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

function run_market_model(data_dir::String, result_dir::String, input_optimizer;
						  return_result::Bool=false, redispatch::Bool=false)
	set_logger()
	@info("Read Model Data..")

	options, data = read_model_data(data_dir)
	data.folders["result_dir"] = result_dir*Dates.format(now(), "dmm_HHMM_SSsss")

	create_folder(data.folders["result_dir"])
	if !("split_timeseries" in keys(options))
		options["split_timeseries"] = false
	end

	if redispatch
		pomato_results = run_market_model_redispatch(data, options, input_optimizer)
	else
		pomato_results = run_market_model(data, options, input_optimizer)
	end

	if return_result
		return pomato_results
	end
	@info("Everything Done!")
end

function run_market_model(data::Data, options::Dict{String, Any}, input_optimizer; 
						  save::Bool=true)

	global optimizer = input_optimizer.Optimizer
	global optimizer_package = input_optimizer
	pomato_results = Dict{String, Result}()
	if options["timeseries"]["split"]
		data_full = deepcopy(data)
		model_horizon_segments = split_timeseries_segments(data_full, options["timeseries"]["market_horizon"])
		for (i, timesteps) in enumerate(model_horizon_segments)
			data = deepcopy(data_full)
			data.t = data.t[timesteps]
			set_model_horizon!(data, i)
			@info("Initializing Market Model for timestep $(data.t[1].name)...")
			pomato_results[data.t[1].name] = market_model(data, options).result
		end
	else
		pomato_results[data.t[1].name] = market_model(data, options).result
	end
	if save
		save_result(concat_results(pomato_results), data.folders["result_dir"])
	end
	return pomato_results
end

function run_market_model_redispatch(data::Data, options::Dict{String, Any}, input_optimizer)

	global optimizer = input_optimizer.Optimizer
	global optimizer_package = input_optimizer

	market_result = concat_results(run_market_model(data, options, input_optimizer, save=false))
	redispatch_results = redispatch_model(market_result, data, options)
	for result in keys(redispatch_results)
		save_result(redispatch_results[result], data.folders["result_dir"]*"_"*result)
	end
	@info("Everything Done!")
	return redispatch_results
end
