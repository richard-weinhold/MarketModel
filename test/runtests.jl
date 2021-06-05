
# cd(pwd()*"/test")
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
			obj_market_result = 1.139441150816381e6
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

		obj_market_result = 1.139441150816381e6
		@test result["market_results"].misc_results["Objective Value"] ≈ obj_market_result rtol=0.01

		for r in keys(result)
			@test result[r].misc_results["Solve Status"] == MarketModel.MOI.OPTIMAL
			@test result[r].misc_results["Objective Value"]*1.01 > obj_market_result
		end
	end
end

@testset "DE" begin
    @testset "Realtime and Dayahead" begin
        data_dir = cd(pwd, "..")*"/examples/de_testing/"
        options, data = MarketModel.read_model_data(data_dir);
        MarketModel.set_da_timeseries!(data)
        MarketModel.set_rt_timeseries!(data)
    end
	@testset "Basic MarketModel" begin
		optimizer_package = Clp
		data_dir = cd(pwd, "..")*"/examples/de_testing/"
		result_dir = cd(pwd, "..")*"/examples/results/"
		result = MarketModel.run_market_model(data_dir, result_dir, optimizer_package, return_result=true)
        @test sum(result["t0001"].H[:, :H]) == 40
	end
end

@testset "Misc Constraints" begin
    @testset "Add PTDF Constraints" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        options, data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data, options);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_dclf_ptdf_constraints!(pomato);
    end
    @testset "Add PTDF Constraints line subset" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        options, data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data, options);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_dclf_ptdf_constraints!(pomato, [1,2]);
    end
    @testset "Add NTC Constraints" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        options, data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data, options);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_ntc_constraints!(pomato);
    end
    @testset "Add NTC Constraints zone subset" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        options, data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data, options);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_ntc_constraints!(pomato, [1]);
    end
    @testset "Add Chance Constraints" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        options, data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data, options);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_chance_constraints!(pomato);
    end
end
