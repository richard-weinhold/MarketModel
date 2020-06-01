
include("../src/MarketModel.jl")
import .MarketModel
using Test, Logging
using Clp

ConsoleLogger(stdout, Logging.Error) |> global_logger#

@testset "All" begin
	@testset "Basic MarketModel" begin
		optimizer = Clp.Optimizer
		data_dir = cd(pwd, "..")*"/examples/nrel_118/"
		result_dir = cd(pwd, "..")*"/examples/results/"
		result = MarketModel.run_market_model(data_dir, result_dir, optimizer, return_result=true)
		for r in keys(result)
			obj_market_result = 993555.8332305112
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
			@test result[r].misc_results["Objective Value"] ≈ obj_market_result rtol=0.01
		end
	end

	@testset "Basic Redispatch" begin
		optimizer = Clp.Optimizer
		data_dir = cd(pwd, "..")*"/examples/nrel_118/"
		result_dir = cd(pwd, "..")*"/examples/results/"
		result = MarketModel.run_market_model(data_dir, result_dir, optimizer,
			redispatch=true, return_result=true)

		obj_market_result = 993555.8332305112
		@test result["market_results"].misc_results["Objective Value"] ≈ obj_market_result rtol=0.01

		for r in keys(result)
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
			@test result[r].misc_results["Objective Value"]*1.01 > obj_market_result
		end
	end
end


#
# ["redispatch_R3"].misc_results["Objective Value"]
#
# data_dir = cd(pwd, "..")*"/examples/nrel_118/"
# options, data = MarketModel.read_model_data(data_dir)
# data.folders["result_dir"] = pwd()*"/examples/results/testtest"
# pomato = MarketModel.POMATO(MarketModel.Model(), data, options)
# MarketModel.add_optimizer!(pomato)
# MarketModel.add_variables_expressions!(pomato)
# pomato.model
# @test MarketModel.JuMP.num_variables(pomato.model) == (211 + 118*3 + 3*3)*24
#
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
