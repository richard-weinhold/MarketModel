
include("../src/MarketModel.jl")
import .MarketModel
using Test, Logging
using Clp

ConsoleLogger(stdout, Logging.Info) |> global_logger

# %%
@testset "All" begin
	@testset "Basic MarketModel" begin
		optimizer_package = Clp
		data_dir = cd(pwd, "..")*"/examples/nrel_118/"
		result_dir = cd(pwd, "..")*"/examples/results/"
		result = MarketModel.run_market_model(data_dir, result_dir, optimizer_package, return_result=true)
		for r in keys(result)
			obj_market_result = 993555.8332305112
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
			@test result[r].misc_results["Objective Value"] ≈ obj_market_result rtol=0.01
		end
	end

	@testset "Basic Redispatch" begin
		optimizer_package = Clp
		data_dir = cd(pwd, "..")*"/examples/nrel_118/"
		result_dir = cd(pwd, "..")*"/examples/results/"
		result = MarketModel.run_market_model(data_dir, result_dir, optimizer_package,
			redispatch=true, return_result=true)

		obj_market_result = 993555.8332305112
		@test result["market_results"].misc_results["Objective Value"] ≈ obj_market_result rtol=0.01

		for r in keys(result)
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
			@test result[r].misc_results["Objective Value"]*1.01 > obj_market_result
		end
	end
end


# %%
# data_dir = "C:/Users/riw/tubCloud/Uni/Market_Tool/pomato_studies/data_temp/julia_files/data/"
# #
# options, data = MarketModel.read_model_data(data_dir)
# options["redispatch"]["zonal_redispatch"] = false
# options["infeasibility"]["electricity"]["bound"] = 0
# optimizer = Gurobi.Optimizer
# data.folders["result_dir"] = cd(pwd, "..")*"/examples/results/"
# # # result = MarketModel.run_market_model(data_dir, result_dir, optimizer,
# # 	# redispatch=true, return_result=true)
# #
# result = MarketModel.run_market_model_redispatch(data, options, optimizer)
#
# maximum(result["redispatch_R1_R2_R3"].INFEAS_EL_N_POS[:, :INFEAS_EL_N_POS])
# maximum(result["redispatch_R1_R2_R3"].INFEAS_EL_N_NEG[:, :INFEAS_EL_N_NEG])
# #
#
# MarketModel.set_global_optimizer(Gurobi.Optimizer)
# pomato = MarketModel.market_model(data, options)
# redispatch_results = Dict{String, MarketModel.Result}()
# redispatch_results["market_results"] = pomato.result
#
# market_result = Dict{String, Array{Float64, 2}}()
# market_result["g_market"] = MarketModel.value.(pomato.model[:G])
# market_result["d_es_market"] = MarketModel.value.(pomato.model[:D_es])
# market_result["d_ph_market"] = MarketModel.value.(pomato.model[:D_ph])
# market_result["infeas_pos_market"] = MarketModel.value.(pomato.model[:INFEAS_EL_N_POS])
# market_result["infeas_neg_market"] = MarketModel.value.(pomato.model[:INFEAS_EL_N_NEG])
#
# MarketModel.load_redispatch_grid!(pomato)
# data_copy = deepcopy(pomato.data)
# market_result_copy = deepcopy(market_result)
#
# if options["redispatch"]["zonal_redispatch"]
# 	redispatch_zones = [[zone] for zone in options["redispatch"]["zones"]]
# else
# 	 redispatch_zones = [convert(Vector{String}, vcat(options["redispatch"]["zones"]))]
# end
#
# # for zones in redispatch_zones
# zones = redispatch_zones[1]
# tmp_results = Dict{String, MarketModel.Result}()
# # for timesteps in [t.index:t.index for t in data_copy.t]
# model_horizon_segments = MarketModel.split_timeseries_segments(data, options["timeseries"]["redispatch_horizon"])
# timesteps = model_horizon_segments[1]
# # data = deepcopy(data_copy)
# data.t = data.t[timesteps]
# MarketModel.set_model_horizon!(data)
# for key in keys(market_result)
# 	market_result[key] = market_result[key][timesteps, :]
# end
# @info("Initializing Redispatch Model for zones $(zones) Timestep $(data.t[1].name)")
#
# options["infeasibility"]["electricity"]["bound"] = 0
# pomato = MarketModel.POMATO(MarketModel.Model(), data, options)
# # MOI.set(pomato.model, MOI.Silent(), false)
# MarketModel.add_optimizer!(pomato);
# MarketModel.redispatch_model!(pomato, market_result, zones);
# MarketModel.add_curtailment_constraints!(pomato, zones);
# MarketModel.add_electricity_energy_balance!(pomato);
#
# @info("Solving...")
# t_start = time_ns()
# MarketModel.JuMP.optimize!(pomato.model)
# @info("Objective: $(MarketModel.JuMP.objective_value(pomato.model))")
# @info("Objective: $(JuMP.termination_status(pomato.model))")
# t_elapsed = time_ns() - t_start
# @info("Solvetime: $(round(t_elapsed*1e-9, digits=2)) seconds")
#


# %%
# obj_market_result = 993555.8332305112
# @test result["market_results"].m isc_results["Objective Value"] ≈ obj_market_result rtol=0.01

# pomato = MarketModel.POMATO(MarketModel.Model(), data, options)


# data.folders["result_dir"] = pwd()*"/examples/results/testtest"
# pomato = MarketModel.POMATO(MarketModel.Model(), data, options)
# MarketModel.add_optimizer!(pomato)
# MarketModel.add_variables_expressions!(pomato)
# pomato.model
# @test MarketModel.JuMP.num_variables(pomato.model) == (211 + 118*3 + 3*3)*24

# include("../src/MarketModel.jl")
# import .MarketModel

# data_dir = cd(pwd, "..")*"/examples/nrel_test/"

# options, data = MarketModel.read_model_data(data_dir)
# options["redispatch"]["horizon"] = 24
# options["redispatch"]["zonal_redispatch"] = false
# optimizer = Clp.Optimizer
#
# data.folders["result_dir"] = pwd()*"/examples/results/testtest"
# result = MarketModel.run_market_model_redispatch(data, options, optimizer)
#
# r1 = result["market_results"]
# # 3.1856171112844176e6
# r2 = result["redispatch_R1_R2_R3"]
# # 3.1830242454394284e6
#
# r1.misc_results["COST_G"] - r2.misc_results["COST_G"]
#
# r1.misc_results["Objective Value"] - r2.misc_results["Objective Value"]

# MarketModel.add_electricity_generation_constraints!(pomato)
# MarketModel.add_electricity_storage_constraints!(pomato)
#
# MarketModel.add_curtailment_constraints!(pomato)
#
# MarketModel.add_ntc_constraints!(pomato)
#
# if in(pomato.options["type"] , ["zonal", "cbco_zonal"])
# 	@info("Adding FlowBased Constraints...")
# 	MarketModel.add_flowbased_constraints!(pomato)
# end
#
# MarketModel.add_dclf_constraints!(pomato)
#
# if (pomato.options["constrain_nex"])
# 	@info("Adding NEX Constraints...")
# 	MarketModel.add_net_position_constraints!(pomato)
# end
#
# MarketModel.add_electricity_energy_balance!(pomato)
#
# @info("Adding Objective Function...")
# MarketModel.add_objective!(pomato)
#
# @time MarketModel.JuMP.optimize!(pomato.model)
# @info("Objective: $(MarketModel.JuMP.objective_value(pomato.model))")
# @info("Objective: $(MarketModel.JuMP.termination_status(pomato.model))")
#
# MarketModel.add_result!(pomato)
# @info("Model Done!")

# "OPTIMAL" in
# for r in t
# 	println(r)
# end
#
# keys(t)
# t["t0001"].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
# result
# options, data = MarketModel.read_model_data(data_dir)
# data.folders["result_dir"] = pwd()*"/examples/results/testtest"
#
# t = MarketModel.run_market_model(data, options, optimizer)
