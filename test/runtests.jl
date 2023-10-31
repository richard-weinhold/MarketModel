
include("../src/MarketModel.jl")
import .MarketModel
using Test, Logging
using Clp
using ECOS

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
        data = MarketModel.read_model_data(data_dir);
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
        data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_dclf_ptdf_constraints!(pomato);
    end
    @testset "Add PTDF Constraints line subset" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_dclf_ptdf_constraints!(pomato, [1,2]);
    end
    @testset "Add NTC Constraints" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_ntc_constraints!(pomato);
    end
    @testset "Add NTC Constraints zone subset" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_ntc_constraints!(pomato, [1]);
    end
    @testset "Add Chance Constraints" begin
        optimizer_package = Clp
        data_dir = cd(pwd, "..")*"/examples/nrel_118/"
        data = MarketModel.read_model_data(data_dir);
        pomato = MarketModel.POMATO(MarketModel.Model(), data);
        MarketModel.add_variables_expressions!(pomato);
        MarketModel.add_chance_constraints!(pomato);
    end
end

@testset "StorageModel" begin
    MarketModel.set_global_optimizer(Clp)	
    data_dir = cd(pwd, "..")*"/examples/nrel_118/"
    data = MarketModel.read_model_data(data_dir);
    data.options["plant_types"]["es"] = ["storage"]
    data.options["timeseries"]["market_horizon"] = 6

    es = MarketModel.Plant(length(data.plants) + 1, "Storage", 1, 5., 0., 0.75, 1., 50., 0., "storage")
    es.storage_capacity = 500
    es.d_max = 50
    es.inflow = zeros(24)
    push!(data.plants, es)

    @test !(isdefined(data.plants[end], :storage_level_start))

    MarketModel.set_storage_levels!(data)

    @test isdefined(data.plants[end], :storage_level_start)
    @test isdefined(data.plants[end], :storage_level_end)
    
    data_full = deepcopy(data)
    model_horizon_segments = MarketModel.split_timeseries_segments(
        data_full, data.options["timeseries"]["market_horizon"]
    )
    @test length(data.plants[end].storage_level_end) == length(model_horizon_segments)
    for (i, timesteps) in enumerate(model_horizon_segments)
        data = deepcopy(data_full)
        data.t = data.t[timesteps]
        MarketModel.set_model_horizon!(data, i)
        @test length(data.t) == 6
    end
end
