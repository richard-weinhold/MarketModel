"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Reading in data as RAW struct. 
"""

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
    storage_level::DataFrame
    reference_flows::DataFrame
    lines::DataFrame
    grid::DataFrame
    redispatch_grid::DataFrame
    contingency_groups::Dict{String, Array{String}}
    slack_zones::DataFrame

    function RAW(data_dir)
        raw = new()
        raw.options = JSON.parsefile(data_dir*"options.json"; dicttype=Dict)
        raw.plant_types = raw.options["plant_types"]

        raw.nodes = DataFrame(CSV.File(data_dir*"nodes.csv"), copycols=true)
        raw.nodes[!, :int_idx] = collect(1:nrow(raw.nodes))
        raw.zones = DataFrame(CSV.File(data_dir*"zones.csv"), copycols=true)
        raw.zones[!, :int_idx] = collect(1:nrow(raw.zones))
        raw.heatareas = DataFrame(CSV.File(data_dir*"heatareas.csv"), copycols=true)
        raw.heatareas[!, :int_idx] = collect(1:nrow(raw.heatareas))
        plants = DataFrame(CSV.File(data_dir*"plants.csv"), copycols=true)
        raw.plants =  filter(row -> !(row[:plant_type] in raw.plant_types["ts"]), plants)
        raw.res_plants =  filter(row -> row[:plant_type] in raw.plant_types["ts"], plants)
        raw.plants[!, :int_idx] = collect(1:nrow(raw.plants))
        raw.res_plants[!, :int_idx] = collect(1:nrow(raw.res_plants))
        raw.availability = DataFrame(CSV.File(data_dir*"availability.csv"), copycols=true)
        raw.demand_el = DataFrame(CSV.File(data_dir*"demand_el.csv"), copycols=true)
        raw.demand_h = DataFrame(CSV.File(data_dir*"demand_h.csv"), copycols=true)
        raw.dc_lines = DataFrame(CSV.File(data_dir*"dclines.csv"), copycols=true)
        raw.ntc = DataFrame(CSV.File(data_dir*"ntc.csv"), copycols=true)
        raw.net_position = DataFrame(CSV.File(data_dir*"net_position.csv"), copycols=true)
        raw.net_export = DataFrame(CSV.File(data_dir*"net_export.csv"), copycols=true)
        
        inflows = DataFrame(CSV.File(data_dir*"inflows.csv"), copycols=true)
        raw.inflows = unstack(inflows, :timestep, :plant, :inflow)

        raw.storage_level = DataFrame(CSV.File(data_dir*"storage_level.csv"), copycols=true)

        raw.lines = DataFrame(CSV.File(data_dir*"lines.csv"), copycols=true)
        raw.grid = DataFrame(CSV.File(data_dir*"grid.csv"), copycols=true)
        raw.redispatch_grid = DataFrame(CSV.File(data_dir*"redispatch_grid.csv"), copycols=true)
        raw.contingency_groups = JSON.parsefile(data_dir*"contingency_groups.json"; dicttype=Dict)
        raw.slack_zones = DataFrame(CSV.File(data_dir*"slack_zones.csv"), copycols=true)
        raw.model_horizon = DataFrame(index=collect(1:size(unique(raw.demand_el[:, :timestep]), 1)),
                                      timesteps=unique(raw.demand_el[:, :timestep]))
        return raw
    end
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
