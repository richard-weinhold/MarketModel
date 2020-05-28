
include("../src/MarketModel.jl")
import .MarketModel
using Test, Logging
using Gurobi, Clp

ConsoleLogger(stdout, Logging.Error) |> global_logger#


@testset "All" begin
	@testset "Basic MarketModel" begin

    optimizer = Clp.Optimizer
    data_dir = pwd()*"/examples/nrel_118/"
    result_dir = pwd()*"/examples/results/"
	end

end

optimizer = Clp.Optimizer
data_dir = pwd()*"/examples/nrel_118/"
result_dir = pwd()*"/examples/results/"
t = MarketModel.run_market_model(data_dir, result_dir, optimizer, return_result=true)

options, data = MarketModel.read_model_data(data_dir)
data.folders["result_dir"] = pwd()*"/examples/results/testtest"

t = MarketModel.run_market_model(data, options, optimizer)
# using Dates
#
#
# options, data = MarketModel.read_model_data(data_dir)
# data.folders = Dict("data_dir" => data_dir,
#                     "result_dir" => result_dir*Dates.format(now(), "dmm_HHMM"))
#
# MarketModel.create_folder(data.folders["result_dir"])
#
#
# # Willy Wonka Manual Adjustments
# # options["infeasibility"]["electricity"]["bound"] = 10000
# # options["type"] = "chance_constrained"
# # data.t = data.t[1:3]
#
# MarketModel.set_model_horizon!(data)
# t = MarketModel.run_market_model(data, options)
#
# MarketModel.Result(t)
#
# result_info = MarketModel.get_result_info(t)
# begin
#     result = MarketModel.Result()
#     for v in keys(result_info)
#         println(v)
#         setfield!(result, v,  MarketModel.model_symbol_to_df(v, result_info, t))
#     end
# end
#
#
# save_result(concat_results(pomato_results), data.folders["result_dir"])
# println("Everything Done!")
# if return_result
#     return pomato_results
# end
