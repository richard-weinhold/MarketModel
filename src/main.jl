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

function run_market_model(data_dir::String, result_dir::String, input_optimizer;
						  return_result::Bool=false, redispatch::Bool=false)
	set_logger()
	@info("Read Model Data..")

	options, data = read_model_data(data_dir)
	data.folders["result_dir"] = result_dir*Dates.format(now(), "dmm_HHMM")

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

function run_market_model(data::Data, options::Dict{String, Any}, input_optimizer)

	global optimizer = input_optimizer
	pomato_results = Dict{String, Result}()
	if options["split_timeseries"]
		data_full = deepcopy(data)
		for timesteps in [t.index:t.index for t in data.t]
			data = deepcopy(data_full)
			data.t = data.t[timesteps]
			set_model_horizon!(data)
			@info("Initializing Market Model for timestep $(data.t[1].name)...")
			pomato_results[data.t[1].name] = run_market_model(data, options).result
		end
	else
		pomato_results[data.t[1].name] = run_market_model(data, options).result
	end
	save_result(concat_results(pomato_results), data.folders["result_dir"])
	return pomato_results
end

function run_market_model_redispatch(data::Data, options::Dict{String, Any}, input_optimizer)

	global optimizer = input_optimizer
	pomato_results = run_redispatch_model(data, options)
	for result in keys(pomato_results)
		save_result(pomato_results[result], data.folders["result_dir"]*"_"*result)
	end
	@info("Everything Done!")
	return pomato_results
end
