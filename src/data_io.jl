
mutable struct RAW
    options::Dict{String, Any}
    model_horizon::DataFrame
    plant_types::Dict{String, Any}
    nodes::DataFrame
    zones::DataFrame
    heatareas::DataFrame
    plants::DataFrame
    res_plants::DataFrame
    availability::DataFrame
    demand_el::DataFrame
    demand_h::DataFrame
    dc_lines::DataFrame
    ntc::DataFrame
    net_position::DataFrame
    net_export::DataFrame
    inflows::DataFrame
    reference_flows::DataFrame
    grid::DataFrame
    slack_zones::DataFrame

    function RAW(data_dir)
        raw = new()
        raw.options = JSON.parsefile(data_dir*"options.json"; dicttype=Dict)
        raw.plant_types = raw.options["plant_types"]

        raw.nodes = DataFrame!(CSV.File(data_dir*"nodes.csv"))
        raw.nodes[!, :int_idx] = collect(1:nrow(raw.nodes))
        raw.zones = DataFrame!(CSV.File(data_dir*"zones.csv"))
        raw.zones[!, :int_idx] = collect(1:nrow(raw.zones))
        raw.heatareas = DataFrame!(CSV.File(data_dir*"heatareas.csv"))
        raw.heatareas[!, :int_idx] = collect(1:nrow(raw.heatareas))
        plants = DataFrame!(CSV.File(data_dir*"plants.csv"))
        raw.plants =  filter(row -> !(row[:plant_type] in raw.plant_types["ts"]), plants)
        raw.res_plants =  filter(row -> row[:plant_type] in raw.plant_types["ts"], plants)
        raw.plants[!, :int_idx] = collect(1:nrow(raw.plants))
        raw.res_plants[!, :int_idx] = collect(1:nrow(raw.res_plants))
        raw.availability = DataFrame!(CSV.File(data_dir*"availability.csv"))
        raw.demand_el = DataFrame!(CSV.File(data_dir*"demand_el.csv"))
        raw.demand_h = DataFrame!(CSV.File(data_dir*"demand_h.csv"))
        raw.dc_lines = DataFrame!(CSV.File(data_dir*"dclines.csv"))
        raw.ntc = DataFrame!(CSV.File(data_dir*"ntc.csv"))
        raw.net_position = DataFrame!(CSV.File(data_dir*"net_position.csv"))
        raw.net_export = DataFrame!(CSV.File(data_dir*"net_export.csv"))
        raw.inflows = DataFrame!(CSV.File(data_dir*"inflows.csv"))
        # raw.reference_flows = DataFrame!(CSV.File(data_dir*"reference_flows.csv");)
        raw.grid = DataFrame!(CSV.File(data_dir*"grid.csv"))
        raw.slack_zones = DataFrame!(CSV.File(data_dir*"slack_zones.csv"))
        raw.model_horizon = DataFrame(index=collect(1:size(unique(raw.demand_el[:, :timestep]), 1)),
                                      timesteps=unique(raw.demand_el[:, :timestep]))
        return raw
    end
end

function load_redispatch_grid!(pomato::POMATO)
    grid = Vector{Grid}()
    grid_data = DataFrame!(CSV.File(pomato.data.folders["data_dir"]*"redispatch_grid.csv"))
    # grid_data = CSV.read(pomato.data.folders["data_dir"]*"redispatch_nodal.csv")
    for cbco in 1:nrow(grid_data)
        index = cbco
        name = grid_data[cbco, :index]
        ptdf = [grid_data[cbco, Symbol(node.name)] for node in pomato.data.nodes]
        ram = grid_data[cbco, :ram]*1.
        newcbco = Grid(index, name, ptdf, ram)
        newcbco.zone_i = coalesce(grid_data[cbco, :zone_i], nothing)
        newcbco.zone_j = coalesce(grid_data[cbco, :zone_j], nothing)
        push!(grid, newcbco)
    end
    pomato.data.grid = grid
    return grid
end

function save_result(result::Result, folder::String)
	create_folder(folder)
	@info("Saving results to folder "*folder*"...")
	for (field, field_type) in zip(fieldnames(Result), fieldtypes(Result))
		if field_type == DataFrame
			CSV.write(folder*"/"*String(field)*".csv",
					  DataFrame(getfield(result, field)))
		elseif field_type == Dict
			open(folder*"/"*String(field)*".json", "w") do f
				# write(f, JSON.json(getfield(result, field), 2))
				JSON.print(f, getfield(result, field), 2)
			end
		else
			@info(field, " not Dict or DataFrame, cant save!")
		end
	end
	@info("All Results Saved!")
end

function create_folder(result_folder)
	if !isdir(result_folder)
	    @info("Creating Results Folder")
		mkpath(result_folder)
	end
end
