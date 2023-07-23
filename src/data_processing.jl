"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Processing of input data and populating data structs.
"""

function read_model_data(data_dir::String)
    @info("Reading Model Data from: $(data_dir)")
    raw = RAW(data_dir)

    nodes = populate_nodes(raw)
    zones = populate_zones(raw, nodes)
    heatareas = populate_heatareas(raw)
    plants = populate_plants(raw)
    res_plants = populate_res_plants(raw)
    dc_lines = populate_dclines(raw)
    
    lines = populate_lines(raw, nodes)
    contingencies, redispatch_contingencies = populate_network(raw, lines, nodes, zones)

    timesteps = populate_timesteps(raw)
    data = Data(nodes, zones, heatareas, plants, res_plants, lines, 
                contingencies, redispatch_contingencies, dc_lines, timesteps)
    data.folders = Dict("data_dir" => data_dir)
    data.options = raw.options
    raw = nothing
    @info("Data Prepared")
    return data
end # end of function

function populate_timesteps(raw::RAW)
    timesteps = Vector{Timestep}()
    for t in 1:nrow(raw.model_horizon)
        index = t
        name = raw.model_horizon[t, :timesteps]
        push!(timesteps, Timestep(index, name))
    end
    return timesteps
end

function populate_zones(raw::RAW, nodes::Vector{Node})
    zones = Vector{Zone}()
    for z in 1:nrow(raw.zones)
        index = z
        name = raw.zones[z, :index]
        nodes_idx = raw.nodes[raw.nodes[:, :zone] .== name, :int_idx]
        nodes_name = raw.nodes[raw.nodes[:, :zone] .== name, :index]
        plants = filter(row -> row[:node] in nodes_name, raw.plants)[:, :int_idx]
        res_plants = filter(row -> row[:node] in nodes_name, raw.res_plants)[:, :int_idx]
        demand = sum(nodes[i].demand for i in nodes_idx)
        # demand = combine(groupby(filter(col -> col[:node] in nodes_name, raw.demand_el), :timestep, sort=true), :demand_el => sum)
        newz = Zone(index, name, demand, nodes_idx, plants, res_plants)
        if (size(raw.ntc, 2) > 1)
            ntc = filter(row -> row[:zone_i] == name, raw.ntc)
            newz.ntc = [zone in ntc[:, :zone_j] ?
                        ntc[ntc[:, :zone_j] .== zone, :ntc][1] :
                        0 for zone in raw.zones[:, :index]]
        end

        newz.net_export = sum(nodes[i].net_export for i in nodes_idx)
        net_position = combine(groupby(filter(col -> col[:zone] == name, raw.net_position), :timestep, sort=true), :net_position => sum)
        if size(net_position, 1) > 0
            newz.net_position = net_position[:, :net_position_sum]
        end
        push!(zones, newz)
    end
    return zones
end

function populate_nodes(raw::RAW)
    nodes = Vector{Node}()
    demand = unstack(raw.demand_el, :timestep, :node, :demand_el)
    net_export = unstack(raw.net_export, :timestep, :node, :net_export)

    for n in 1:nrow(raw.nodes)
        index = n
        name = raw.nodes[n, :index]
        slack = raw.nodes[n, :slack]
        zone_name = raw.nodes[n, :zone]
        zone_idx = raw.zones[raw.zones[:, :index] .== zone_name, :int_idx][1]
        plants = raw.plants[raw.plants[:, :node] .== name, :int_idx]
        res_plants = raw.res_plants[raw.res_plants[:, :node] .== name, :int_idx]
        newn = Node(index, name, zone_idx, demand[:, name], slack, plants, res_plants)
        if slack
            # newn.slack_zone = slack_zones[index]
            slack_zone = raw.slack_zones[:, :index][raw.slack_zones[:, Symbol(name)] .== 1]
            newn.slack_zone = filter(col -> col[:index] in slack_zone, raw.nodes)[:, :int_idx]
        end
        newn.net_export = (name in names(net_export) ? net_export[:, name] :
                           zeros(length(raw.model_horizon[:, :timesteps])))

        if "demand_el_da" in string.(names(raw.demand_el))
            demand_el_da = combine(groupby(filter(col -> col[:node] == name, raw.demand_el),
                                      :timestep, sort=true), :demand_el_da => sum)
            newn.demand_da =  demand_el_da[:, :demand_el_da_sum]
        end
        push!(nodes, newn)
    end
    return nodes
end

function populate_heatareas(raw::RAW)
    heatareas = Vector{Heatarea}()
    for h in 1:nrow(raw.heatareas)
        index = h
        name = raw.heatareas[h, :index]
        demand = combine(groupby(filter(col -> col[:heatarea] == name, raw.demand_h), :timestep, sort=true), :demand_h => sum)
        plants = raw.plants[(raw.plants[:, :heatarea] .=== name).&(raw.plants[:, :h_max] .> 0), :int_idx]
        res_plants = raw.res_plants[(raw.res_plants[:, :heatarea] .=== name).&(raw.res_plants[:, :h_max] .> 0), :int_idx]
        newh = Heatarea(index, name, demand[:, :demand_h_sum], plants, res_plants)
        push!(heatareas, newh)
    end
    return heatareas
end

function populate_plants(raw::RAW)
    plants =  Vector{Plant}()
    es_plants = union(raw.plant_types["hs"], raw.plant_types["es"])
    for p in 1:nrow(raw.plants)
        index = p
        name = string(raw.plants[p, :index])
        node_name = raw.plants[p, :node]
        node_idx = raw.nodes[raw.nodes[:, :index] .== node_name, :int_idx][1]
        eta = raw.plants[p, :eta]*1.
        availability = ("availability" in names(raw.plants) ? raw.plants[p, :availability]*1. : 1.)
        g_max = raw.plants[p, :g_max]*1.
        h_max = raw.plants[p, :h_max]*1.
        mc_el = raw.plants[p, :mc_el]*1.
        mc_heat = raw.plants[p, :mc_heat]*1.
        plant_type = raw.plants[p, :plant_type]
        newp = Plant(index, name, node_idx, mc_el,
                     mc_heat, eta, availability, g_max, h_max, plant_type)

        if plant_type in union(raw.plant_types["hs"], raw.plant_types["es"])
            newp.inflow = (
                name in names(raw.inflows) ? 
                raw.inflows[:,Symbol(name)] :
                zeros(length(raw.model_horizon[:, :timesteps]))
            )
            newp.storage_level_start = raw.storage_level[raw.storage_level[:, :plant] .== name, :storage_start]
            newp.storage_level_end = raw.storage_level[raw.storage_level[:, :plant] .== name, :storage_end]
            newp.storage_start = raw.options["storages"]["storage_start"]
            newp.storage_end = raw.options["storages"]["storage_end"]
            newp.storage_capacity = raw.plants[p, :storage_capacity]
            newp.d_max = raw.plants[p, :d_max]
        end
        push!(plants, newp)
    end
    return plants
end

function populate_res_plants(raw::RAW)
    res_plants = Vector{Renewables}()
    if size(raw.availability, 1) > 0
        availability = unstack(raw.availability, :timestep, :plant, :availability)
    end
    sigma_factor = raw.options["chance_constrained"]["percent_std"]
    for res in 1:nrow(raw.res_plants)
        index = res
        name = string(raw.res_plants[res, :index])
        node_name = raw.res_plants[res, :node]
        node_idx = raw.nodes[raw.nodes[:, :index] .== node_name, :int_idx][1]
        g_max = raw.res_plants[res, :g_max]*1.
        h_max = raw.res_plants[res, :h_max]*1.
        mc_el = raw.res_plants[res, :mc_el]*1.
        mc_heat = raw.res_plants[res, :mc_heat]*1.
        plant_type = raw.res_plants[res, :plant_type]
        newres = Renewables(
            index, name, g_max, h_max, mc_el, mc_heat, availability[:, name], node_idx, plant_type, sigma_factor
        )

        if "availability_da" in string.(names(raw.availability))
            availability_da = filter(col -> col[:plant] == name, raw.availability)[:, :availability_da]
            
            newres.mu_da = availability_da * g_max
            newres.sigma_da = newres.mu_da * newres.sigma_factor 
            newres.mu_heat_da = availability_da * h_max
            newres.sigma_heat_da = newres.mu_heat_da * newres.sigma_factor
        end
        push!(res_plants, newres)
    end
    return res_plants
end

function populate_dclines(raw::RAW)
    dc_lines = Vector{DC_Line}()
    for dc in 1:nrow(raw.dc_lines)
        index = dc
        name = raw.dc_lines[dc, :index]
        node_i = raw.dc_lines[dc, :node_i]
        node_j = raw.dc_lines[dc, :node_j]
        node_i_idx = raw.nodes[raw.nodes[:, :index] .== node_i, :int_idx][1]
        node_j_idx = raw.nodes[raw.nodes[:, :index] .== node_j, :int_idx][1]
        capacity = raw.dc_lines[dc, :capacity]*1.
        newdc = DC_Line(index, name, node_i_idx, node_j_idx, capacity)
        push!(dc_lines, newdc)
    end
    return dc_lines
end

function populate_lines(raw::RAW, nodes::Vector{Node})
    lines = Vector{Line}()
    for line in 1:nrow(raw.lines)
        index = line
        name = raw.lines[line, :index]
        capacity = 1. * raw.lines[line, :capacity]
        x = 1. * raw.lines[line, :x_pu]
        b = 1. / raw.lines[line, :x_pu]
        incidence = zeros(Int, length(nodes))
        incidence[[findfirst(n -> n.name == raw.lines[line, :node_i], nodes), 
                findfirst(n -> n.name == raw.lines[line, :node_j], nodes)]] = [1,-1]

        zone_i = nodes[findfirst(n -> n.name == raw.lines[line, :node_i], nodes)].zone
        zone_j = nodes[findfirst(n -> n.name == raw.lines[line, :node_j], nodes)].zone
        newline = Line(index, name, b, capacity, zone_i, zone_j, incidence)
        push!(lines, newline)
    end
    return lines
end

function populate_contingencies(grid_data::DataFrame, raw::RAW, lines::Vector{Line}, nodes::Vector{Node})
    contingencies = Vector{Contingency}()
    for contingency in unique(grid_data[:, :co])
        name = contingency
        cb = findall(l -> l.name in grid_data[grid_data[:, :co] .== contingency, :cb], lines)
        co = (contingency in keys(raw.contingency_groups) ? findall(l -> l.name in raw.contingency_groups[contingency], lines) : Vector{Int}())
        ptdf = Array(grid_data[grid_data[:, :co] .== contingency, [Symbol(n.name) for n in nodes]])
        ram = Vector(grid_data[grid_data[:, :co] .== contingency, :ram]).*1.
        newcontingency = Contingency(name, cb, co, ptdf, ram)
        push!(contingencies, newcontingency)
    end
    return contingencies
end

function populate_contingencies(grid_data::DataFrame, raw::RAW, lines::Vector{Line}, zones::Vector{Zone})
    contingencies = Vector{Contingency}()
    for contingency in unique(grid_data[:, :co])
        name = contingency
        tmp_co = grid_data[grid_data[:, :co] .== contingency, :cb]
        cb = [findfirst(l -> l.name==line, lines) for line in tmp_co]
        co = (contingency in keys(raw.contingency_groups) ? findall(l -> l.name in raw.contingency_groups[contingency], lines) : Vector{Int}())
        ptdf = Array(grid_data[grid_data[:, :co] .== contingency, [Symbol(z.name) for z in zones]])
        ram = Vector(grid_data[grid_data[:, :co] .== contingency, :ram]).*1.
        newcontingency = Contingency(name, cb, co, ptdf, ram)
        push!(contingencies, newcontingency)
    end
    return contingencies
end

function populate_network(raw::RAW, lines::Vector{Line}, nodes::Vector{Node}, zones::Vector{Zone})

    if raw.options["type"] in ["fbmc"]
        if "timestep" in names(raw.grid)
            tmp_contingencies = []
            for timestep in raw.model_horizon[:, :timesteps]
                tmp = populate_contingencies(raw.grid[raw.grid[:, :timestep] .== timestep, :],
                                             raw, lines, zones)
                map(c -> c.timestep = timestep, tmp)
                push!(tmp_contingencies, tmp)
            end
            contingencies = vcat(tmp_contingencies...)
        else
            contingencies = populate_contingencies(raw.grid, raw, lines, zones)
        end
    else # raw.options["type"] in ["nodal", "cbco_nodal"]
        contingencies = populate_contingencies(raw.grid, raw, lines, nodes)
    end
    redispatch_contingencies = populate_contingencies(raw.redispatch_grid, raw, lines, nodes)
    return contingencies, redispatch_contingencies
end

function set_model_horizon!(data::Data, split::Int)
	timesteps = [t.index for t in data.t]
	for n in data.nodes
		n.demand = n.demand[timesteps]
		n.net_export = n.net_export[timesteps]
	end
	for z in data.zones
		z.demand = z.demand[timesteps]
		if any([isdefined(z, :net_position) for z in data.zones])
			z.net_position = z.net_position[timesteps]
		end
		z.net_export = z.net_export[timesteps]
	end
    for p in filter(plant -> isdefined(plant, :storage_level_start), data.plants)
        p.storage_start = p.storage_level_start[split]
        p.storage_end = p.storage_level_end[split]
        p.inflow = p.inflow[timesteps]
    end
    for res in data.renewables
		res.mu = res.mu[timesteps]
		res.mu_heat = res.mu_heat[timesteps]
		res.sigma = res.sigma[timesteps]
		res.sigma_heat = res.sigma_heat[timesteps]
	end
end

function set_model_horizon!(data::Data)
	timesteps = [t.index for t in data.t]
	for n in data.nodes
		n.demand = n.demand[timesteps]
		n.net_export = n.net_export[timesteps]
	end
	for z in data.zones
		z.demand = z.demand[timesteps]
		if any([isdefined(z, :net_position) for z in data.zones])
			z.net_position = z.net_position[timesteps]
		end
		z.net_export = z.net_export[timesteps]
	end
    for res in data.renewables
		res.mu = res.mu[timesteps]
		res.mu_rt = res.mu_rt[timesteps]
        if isdefined(res, :mu_da)
		    res.mu_da = res.mu_da[timesteps]
        end
		res.mu_heat = res.mu_heat[timesteps]
		res.sigma = res.sigma[timesteps]
		res.sigma_heat = res.sigma_heat[timesteps]
	end
end

function set_da_timeseries!(data::Data)

    if all(isdefined(res, :mu_da) for res in data.renewables)
        @info("Using Day Ahead timeseries for availability")
    	for res in data.renewables
    		res.mu = copy(res.mu_da)
    		res.mu_heat = copy(res.mu_heat_da)
    		res.sigma = copy(res.sigma_da)
    		res.sigma_heat = copy(res.sigma_heat_da)
            res.real_time = false
    	end
    else
        @info("No DA timeseries for availability, using RT data!")
    end
    if all(isdefined(node, :demand_da) for node in data.nodes)
        @info("Using Day Ahead timeseries for demand")
    	for node in data.nodes
    		node.demand = copy(node.demand_da)
            node.real_time = false
    	end
    	for zone in data.zones
    		zone.demand = sum(n.demand for n in data.nodes[zone.nodes])
    	end
    else
        @info("No DA timeseries for demand, using RT data!")
    end
end

function set_rt_timeseries!(data::Data)
    @info("Using Real Time timeseries")
	for res in data.renewables
		res.mu = copy(res.mu_rt)
		res.mu_heat = copy(res.mu_heat_rt)
		res.sigma = copy(res.sigma_rt)
		res.sigma_heat = copy(res.sigma_heat_rt)
        res.real_time = false
	end
	for node in data.nodes
		node.demand = copy(node.demand_rt)
        node.real_time = false
	end
	for zone in data.zones
		zone.demand = sum(n.demand for n in data.nodes[zone.nodes])
	end
end
