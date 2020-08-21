
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
# cd("test")
# pwd()
# include("../src/MarketModel.jl")
# import .MarketModel
# using Test, Logging
# using Clp, Gurobi, Mosek, MosekTools
# using DataFrames

# data_dir = "C:/Users/riw/tubCloud/Uni/Market_Tool/pomato_studies/data_temp/julia_files/data/"
# # options, data = MarketModel.read_model_data(data_dir)
# # raw = MarketModel.RAW(data_dir)

# result_dir = cd(pwd, "..")*"/examples/results/"
# @time result = MarketModel.run_market_model(data_dir, result_dir, Mosek, return_result=true, redispatch=true)

# println("Done")
