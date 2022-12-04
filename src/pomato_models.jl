"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Definition of the models defined by input data and user options. 
"""

function add_optimizer!(pomato::POMATO)
	global optimizer
	global optimizer_package
	set_optimizer(pomato.model, optimizer)
	if string(optimizer_package) == "Gurobi"
		set_optimizer_attributes(
			pomato.model, 
			"LogFile" => pomato.data.folders["result_dir"]*"/log.txt"
		)
		if ("solver_options" in keys(pomato.options["solver"])) && (length(pomato.options["solver"]["solver_options"]) > 0)
			@info("Adding user solver options: ")
			for option in keys(pomato.options["solver"]["solver_options"])
				@info("$(option): $(pomato.options["solver"]["solver_options"][option])")
				set_optimizer_attribute(
					pomato.model, option, pomato.options["solver"]["solver_options"][option])
			end
		else
			@info("Adding default solver options for Gurobi: Method 1, Threads: $(Threads.nthreads() - 2)")
			set_optimizer_attribute(pomato.model, "Method", 1)
			set_optimizer_attribute(pomato.model, "Threads", Threads.nthreads() - 2)
		end
	end
end

function market_model(data)

	pomato = POMATO(Model(), data)

	if pomato.options["timeseries"]["type"] == "da"
		set_da_timeseries!(data)
	end

	add_optimizer!(pomato)

	@info("Adding Variables and Expressions...")
	add_variables_expressions!(pomato)

	@info("Adding Base Model...")
	add_electricity_generation_constraints!(pomato)
	add_electricity_storage_constraints!(pomato)

	if pomato.options["heat"]["include"]
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

	if in(pomato.options["type"] , ["fbmc"])
		if pomato.options["chance_constrained"]["include"]
			@info("Adding FB-Chance Constraints...")
			@time add_chance_constrained_flowbased_constraints!(pomato)
		else
			@info("Adding FlowBased Constraints...")
			add_flowbased_constraints!(pomato)
		end
		non_fb_region = findall(zone -> !(zone.name in pomato.options["fbmc"]["flowbased_region"]), data.zones)
		if length(non_fb_region) > 0
			add_ntc_constraints!(pomato, non_fb_region)
		end
	end

	if in(pomato.options["type"] , ["opf", "scopf"]) 
		if pomato.options["chance_constrained"]["include"]
			@info("Adding Chance Constraints...")
			add_chance_constraints!(pomato)
		else	
			if length(data.contingencies) > 1
				@info("Adding power flow constraints using the PTDF formulation.")
				add_dclf_ptdf_constraints!(pomato)
			else
				@info("Adding power flow constraints using the angle formulation.")
				add_dclf_angle_constraints!(pomato)
				# add_ntc_constraints!(pomato)
				# add_dclf_ptdf_constraints!(pomato)
			end
		end
	end

	add_electricity_energy_balance!(pomato::POMATO)

	@info("Adding Objective Function...")
	add_objective!(pomato);

	@info("Solving...")
	t_start = time_ns()
	@time JuMP.optimize!(pomato.model)
	@info("Termination Status: $(JuMP.termination_status(pomato.model))")
	if JuMP.termination_status(pomato.model) == MOI.INFEASIBLE
		# check_infeasibility(pomato)
		throw("Model is Infeasible. See stored information from Gurobi constraint conflics.")
	elseif JuMP.termination_status(pomato.model) != MOI.OPTIMAL
		@info("Termination Status not optimal, check solution for feasibility.")
	end
	@info("Objective: $(JuMP.objective_value(pomato.model))")
	@info("Optimization Status: $(JuMP.termination_status(pomato.model))")
	t_elapsed = time_ns() - t_start
	@info("Solvetime: $(round(t_elapsed*1e-9, digits=2)) seconds")
	add_result!(pomato)
	@info("Model Done!")
	return pomato
end

function redispatch_model(market_result::Result, data::Data)

	set_rt_timeseries!(data)
	redispatch_results = Dict{String, Result}()
	redispatch_results["market_results"] = market_result

	market_result_variables = Dict{String, Array{Float64, 2}}()
	mapping_he = findall(plant -> plant.h_max > 0, data.plants)
	es = findall(plant -> plant.plant_type in data.options["plant_types"]["es"], data.plants)
	ph = findall(plant -> plant.plant_type in data.options["plant_types"]["ph"], data.plants[mapping_he])

	market_result_variables["g_market"] = Array(sort(unstack(market_result.G, :t, :p, :G))[:, [p.name for p in data.plants]])
	market_result_variables["curt_market"] = (
		size(market_result.CURT, 1) > 0 ? 
		Array(sort(unstack(market_result.CURT, :t, :p, :CURT))[:, [res.name for res in data.renewables]]) : 
		zeros(length(data.t), length(data.renewables))
	)
	market_result_variables["d_es_market"] = (
		size(market_result.D_es, 1) > 0 ? 
		Array(sort(unstack(market_result.D_es, :t, :p, :D_es))[:, [p.name for p in data.plants[es]]]) : 
		Array{Float64}(undef, length(data.t), 0)
	)
	market_result_variables["d_ph_market"] = (size(market_result.D_ph, 1) > 0 ? 
		Array(sort(unstack(market_result.D_ph, :t, :p, :D_ph))[:, [p.name for p in data.plants[ph]]]) : 
		Array{Float64}(undef, length(data.t), 0)
	)
	market_result_variables["infeas_pos_market"] = Array(sort(unstack(market_result.INFEASIBILITY_EL_POS, :t, :n, :INFEASIBILITY_EL_POS))[:, [n.name for n in data.nodes]])
	market_result_variables["infeas_neg_market"] = Array(sort(unstack(market_result.INFEASIBILITY_EL_NEG, :t, :n, :INFEASIBILITY_EL_NEG))[:, [n.name for n in data.nodes]])

	data.contingencies = data.redispatch_contingencies
	pomato = POMATO(Model(), data)
	
	data_copy = deepcopy(pomato.data)
	market_result_variables_copy = deepcopy(market_result_variables)

	if data.options["redispatch"]["zonal_redispatch"]
		redispatch_zones = [[zone] for zone in data.options["redispatch"]["zones"]]
	else
		 redispatch_zones = [convert(Vector{String}, vcat(data.options["redispatch"]["zones"]))]
	end
	for zones in redispatch_zones
		tmp_results = Dict{String, Result}()
		# for timesteps in [t.index:t.index for t in data_copy.t]
		model_horizon_segments = split_timeseries_segments(
			data, data.options["timeseries"]["redispatch_horizon"]
		)
		for timesteps in model_horizon_segments #[1:1]
			# data = deepcopy(data_copy)
			pomato = solve_redispatch_model(
				deepcopy(data_copy), deepcopy(market_result_variables_copy), timesteps, zones
			)
			
			tmp_results[data_copy.t[timesteps][1].name] = add_result!(pomato)
		end
		redispatch_results["redispatch_"*join(zones, "_")] = concat_results(tmp_results)
	end
	return redispatch_results
end

function solve_redispatch_model(
	data::Data, market_result_variables::Dict{String, Array{Float64, 2}},
	timesteps::UnitRange, redispatch_zones::Vector{String}
	)

	data.t = data.t[timesteps]
	set_model_horizon!(data)
	for key in keys(market_result_variables)
		market_result_variables[key] = market_result_variables[key][timesteps, :]
	end
	@info("Initializing Redispatch Model for zones $(redispatch_zones) Timestep $(data.t[1].name)")

	pomato = POMATO(Model(), data)
	MOI.set(pomato.model, MOI.Silent(), false)
	add_optimizer!(pomato);
	redispatch_model!(pomato, market_result_variables, redispatch_zones);
	add_curtailment_constraints!(pomato, redispatch_zones, market_result_variables["curt_market"]);
	add_electricity_energy_balance!(pomato);
	@info("Solving...")
	t_start = time_ns()
	JuMP.optimize!(pomato.model)
	@info("Objective: $(JuMP.objective_value(pomato.model))")
	@info("Objective: $(JuMP.termination_status(pomato.model))")
	t_elapsed = time_ns() - t_start
	@info("Solvetime: $(round(t_elapsed*1e-9, digits=2)) seconds")

	if JuMP.termination_status(pomato.model) != MOI.OPTIMAL
		# check_infeasibility(pomato)
	end
	return pomato
end
