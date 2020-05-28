
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
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
		end
	end
end

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
