"""
asd

"""

function add_optimizer!(pomato::POMATO)
	global optimizer
	global optimizer_package
	set_optimizer(pomato.model, optimizer)
	if string(optimizer_package) == "Gurobi"
		set_optimizer_attributes(pomato.model, "Method" => 1,
								 "Threads" => Threads.nthreads() - 2,
								 "LogFile" => pomato.data.folders["result_dir"]*"/log.txt")
	end
end

function market_model(data::Data, options::Dict{String, Any})

	pomato = POMATO(Model(), data, options)

	add_optimizer!(pomato)

	@info("Adding Variables and Expressions..")
	add_variables_expressions!(pomato)

	@info("Add Base Model")
	add_electricity_generation_constraints!(pomato)
	add_electricity_storage_constraints!(pomato)

	if pomato.options["heat_model"]
		@info("Adding Heat Model...")
		add_heat_generation_constraints!(pomato)
	end
	if pomato.options["curtailment"]["include"]
		@info("Adding Curtailment...")
		add_curtailment_constraints!(pomato)
	end

	if in(pomato.options["type"] , ["ntc"])
		@info("Adding NTC Constraints...")
		add_ntc_constraints!(pomato)
	end

	if in(pomato.options["type"] , ["zonal", "cbco_zonal"])
		@info("Adding FlowBased Constraints...")
		add_flowbased_constraints!(pomato)
	end

	if in(pomato.options["type"] , ["cbco_nodal", "nodal"]) & !options["chance_constrained"]["include"]
		@info("Adding Load Flow Constraints...")
		add_dclf_constraints!(pomato)
	end

	if options["chance_constrained"]["include"]
		@info("Adding Chance Constraints...")
		@time add_chance_constraints!(pomato)
	end

	if (pomato.options["constrain_nex"])
		@info("Adding NEX Constraints...")
		add_net_position_constraints!(pomato)
	end

	add_electricity_energy_balance!(pomato::POMATO)

	@info("Adding Objective Function...")
	add_objective!(pomato)

	@info("Solving...")
	t_start = time_ns()
	@time JuMP.optimize!(pomato.model)
	@info("Objective: $(JuMP.objective_value(pomato.model))")
	@info("Objective: $(JuMP.termination_status(pomato.model))")
	t_elapsed = time_ns() - t_start
	@info("Solvetime: $(round(t_elapsed*1e-9, digits=2)) seconds")
	if JuMP.termination_status(pomato.model) != MOI.OPTIMAL
		check_infeasibility(pomato.model)
	end
	add_result!(pomato)
	@info("Model Done!")
	return pomato
end

function redispatch_model(data::Data, options::Dict{String, Any})

	set_rt_timeseries!(data)
	pomato = market_model(data, options)
	set_rt_timeseries!(pomato.data)

	redispatch_results = Dict{String, Result}()
	redispatch_results["market_results"] = pomato.result

	market_result = Dict{String, Array{Float64, 2}}()
	market_result["g_market"] = value.(pomato.model[:G])
	market_result["d_es_market"] = value.(pomato.model[:D_es])
	market_result["d_ph_market"] = value.(pomato.model[:D_ph])
	market_result["infeas_pos_market"] = value.(pomato.model[:INFEAS_EL_N_POS])
	market_result["infeas_neg_market"] = value.(pomato.model[:INFEAS_EL_N_NEG])

	load_redispatch_grid!(pomato)
	data_copy = deepcopy(pomato.data)
	market_result_copy = deepcopy(market_result)

	if options["redispatch"]["zonal_redispatch"]
		redispatch_zones = [[zone] for zone in options["redispatch"]["zones"]]
	else
		 redispatch_zones = [convert(Vector{String}, vcat(options["redispatch"]["zones"]))]
	end

	for zones in redispatch_zones
		tmp_results = Dict{String, Result}()
		# for timesteps in [t.index:t.index for t in data_copy.t]
		model_horizon_segments = split_timeseries_segments(data, options["timeseries"]["redispatch_horizon"])
		for timesteps in model_horizon_segments
			# data = deepcopy(data_copy)
			pomato = solve_redispatch_model(deepcopy(data_copy), deepcopy(market_result_copy), options, timesteps, zones)
			tmp_results[data_copy.t[timesteps][1].name] = add_result!(pomato)
		end
		redispatch_results["redispatch_"*join(zones, "_")] = concat_results(tmp_results)
	end
	return redispatch_results
end

function solve_redispatch_model(data::Data, market_result::Dict{String, Array{Float64, 2}},
								options::Dict{String, Any}, timesteps::UnitRange, redispatch_zones::Vector{String})
	data.t = data.t[timesteps]
	set_model_horizon!(data)
	for key in keys(market_result)
		market_result[key] = market_result[key][timesteps, :]
	end
	@info("Initializing Redispatch Model for zones $(redispatch_zones) Timestep $(data.t[1].name)")

	pomato = POMATO(Model(), data, options)
	MOI.set(pomato.model, MOI.Silent(), false)
	add_optimizer!(pomato);
	redispatch_model!(pomato, market_result, redispatch_zones);
	add_curtailment_constraints!(pomato, redispatch_zones);
	add_electricity_energy_balance!(pomato);

	# add_objective!(pomato)
	@info("Solving...")
	t_start = time_ns()
	JuMP.optimize!(pomato.model)
	@info("Objective: $(JuMP.objective_value(pomato.model))")
	@info("Objective: $(JuMP.termination_status(pomato.model))")
	t_elapsed = time_ns() - t_start
	@info("Solvetime: $(round(t_elapsed*1e-9, digits=2)) seconds")

	if JuMP.termination_status(pomato.model) != MOI.OPTIMAL
		check_infeasibility(pomato.model)
	end
	return pomato
end
