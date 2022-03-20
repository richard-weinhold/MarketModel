

mutable struct ZonalGenerationCapacity
    capacity::Vector{Float64}
    mc::Vector{Float64}
    res_infeed::Vector{Float64}
    function ZonalGenerationCapacity(
        plants::Vector{MarketModel.Plant}, 
        res::Vector{MarketModel.Renewables};
        clusters=5)
        zonal_capacity = new()
        if length(plants) > clusters
            R = Clustering.kmeans(reshape([p.mc_el for p in plants], 1, length(plants)), clusters)
            zonal_capacity.capacity = [sum(p -> p.g_max, plants[assignments(R) .== i]) for i in 1:clusters]
            zonal_capacity.mc = R.centers[1,:]
        elseif (0 < length(plants) <= clusters)
            zonal_capacity.capacity = vcat([p.g_max for p in plants], zeros(clusters - length(plants)))
            zonal_capacity.mc = vcat([p.mc_el for p in plants], zeros(clusters - length(plants)))
        else
            zonal_capacity.capacity = zeros(clusters)
            zonal_capacity.mc = zeros(clusters)
        end
        zonal_capacity.res_infeed = sum(res -> res.mu, res)
        return zonal_capacity
    end
end

function solve_storage_model(data::Data)
    # Setup 
    nth_hour = 5
    T = collect(1:nth_hour:length(data.t)-nth_hour)
    Z = 1:length(data.zones)
    C = 1:5
    zonal_capacity = Vector{ZonalGenerationCapacity}()
    for z in Z
        push!(zonal_capacity, ZonalGenerationCapacity(
            filter(p -> p.node in data.zones[z].nodes && !(p.plant_type in data.options["plant_types"]["es"]), data.plants),
            filter(p -> (p.node in data.zones[z].nodes), data.renewables),
        ))
    end
    storages = filter(plant -> plant.plant_type in data.options["plant_types"]["es"], data.plants)
    ES = 1:length(storages)
    small_ralative_capacity = [es.storage_capacity < 12*es.g_max for es in storages]
    small_absolute_capacity = [es.g_max < 25 for es in storages]
    short_term_es = ES[small_ralative_capacity .| small_absolute_capacity]
    T_short_term_es = T[1:floor(Int, 168*1/nth_hour):length(T)]
    @info("Solving Storage Model for $(length(T)) timesteps")
    # Define and solve model
    StorageModel = Model()
    @variable(StorageModel, 0 <= G[t=T, z=Z, c=C] <= zonal_capacity[z].capacity[c]);
    @variable(StorageModel, 0 <= L[t=T, es=ES] <= storages[es].storage_capacity);
    
    @variable(StorageModel, 0 <= G_es[t=T, es=ES] <= storages[es].g_max);
    @variable(StorageModel, 0 <= D_es[t=T, es=ES] <= storages[es].d_max);
    @variable(StorageModel, 0 <= CURT[t=T, z=Z] <= mean(zonal_capacity[z].res_infeed[t:(t+nth_hour-1)]));
    @variable(StorageModel, 0 <= EX[t=T, z=Z, zz=Z] <= data.zones[z].ntc[zz]);
    @variable(StorageModel, 0 <= INF_POS[t=T, z=Z] <= 1e9);
    @variable(StorageModel, 0 <= INF_NEG[t=T, z=Z] <= 1e9);
    @variable(StorageModel, 
        0 <= Dump_Water[t=T, es=ES] <= sum(storages[es].inflow[t:(t+nth_hour-1)])
        # 0 <= Dump_Water[t=T, es=ES]
    );

    @expression(StorageModel, G_es_zone[t=T, z=Z], 
        sum(G_es[t, es] for es in findall(es -> es.node in data.zones[z].nodes, storages))
    );

    @expression(StorageModel, D_es_zone[t=T, z=Z], 
        sum(D_es[t, es] for es in findall(es -> es.node in data.zones[z].nodes, storages))
    );
    
    @expression(StorageModel, COST_G[t=T], 
        sum(G[t,z,c]*zonal_capacity[z].mc[c] for z in Z, c in C)        
    );
    
    @objective(StorageModel, Min, 
        + sum(COST_G)
        + sum(CURT)*data.options["curtailment"]["cost"] + sum(EX)*1.
        + sum(INF_POS)*1e4 + sum(INF_NEG)*1e4
    );
    
    @constraint(StorageModel, EnergyBalance[t=T, z=Z],
        mean(data.zones[z].demand[t:(t+nth_hour-1)]) ==
        sum(G[t,z,c] for c in C) + G_es_zone[t,z] 
        + mean(zonal_capacity[z].res_infeed[t:(t+nth_hour-1)])
        - D_es_zone[t,z]
        - CURT[t,z]
        - sum(EX[t, z, zz] for zz in Z)
        + sum(EX[t, zz, z] for zz in Z)
        + INF_POS[t, z] - INF_NEG[t, z]
    );
    
    storage_start(es) = storages[es].storage_capacity*0.5
    @constraint(StorageModel, [t=T, es=ES],
        L[t, es] == (t>1 ? L[t-nth_hour, es] : storage_start(es)) 
        + sum(storages[es].inflow[t:(t+nth_hour-1)]) 
        - G_es[t, es]*nth_hour - 
        Dump_Water[t, es]*nth_hour + storages[es].eta*D_es[t, es]*nth_hour
    );
    # lower bound on storage level in last timestep
    @constraint(StorageModel, StorageEnd[es=ES],
        L[T[end], es] >= storage_start(es)
    );
    # Short Term storages balance every 4 weeks
    @constraint(StorageModel, ShortTermES[es=short_term_es, t=T_short_term_es],
        L[t, es] == storage_start(es)
    );
    global optimizer
    global optimizer_package
    set_optimizer(StorageModel, optimizer)   
    if string(optimizer_package) == "Gurobi"
        set_optimizer_attribute(StorageModel, "Method", 2)
        set_optimizer_attribute(StorageModel, "Crossover", 0)
    end
    optimize!(StorageModel)
    return StorageModel
end

function set_storage_levels!(data::Data)
    storages = filter(p -> p.plant_type in data.options["plant_types"]["es"], data.plants)
    storage_model = solve_storage_model(data)
    model_horizon_segments = split_timeseries_segments(
        data, data.options["timeseries"]["market_horizon"]
    )
    default_value = data.options["parameters"]["storage_start"]
    default_level = [default_value for i in 1:length(model_horizon_segments)]
    
    ES = 1:length(storages)
    storage_model_level = value.(storage_model[:L]).data[:, ES]
    storage_capacity = [storages[es].storage_capacity for es in ES]'
    relative_storage_level =  storage_model_level ./ storage_capacity
    cycles = sum(abs.(diff(relative_storage_level, dims=1)), dims=1)[1, :]
    
    # Step 1) Clasifucation: Short Term Storages 
    # fully charge and discharge the storage more than 8 times a year.
    many_cycles = (cycles .> 8) 
    small_ralative_capacity = [es.storage_capacity < 12*es.g_max for es in storages]
    small_absolute_capacity = [es.g_max < 25 for es in storages]
    for es in storages[(many_cycles .| small_ralative_capacity .| small_absolute_capacity)]
        es.storage_level_start = default_level
        es.storage_level_end = default_level
    end 
    
    # Seasonal levels for all other storages
    T = storage_model[:L].axes[1]
    window = 7
    for es in ES[(.!many_cycles .& .!small_ralative_capacity .& .!small_absolute_capacity)]
        level_rolling = RollingFunctions.rollmean(relative_storage_level[:,es], window)
        inter = Interpolations.LinearInterpolation(T[window:end], level_rolling)
        storage_level_start = [inter(s.start) for s in model_horizon_segments[2:end]]
        pushfirst!(storage_level_start, 0.5)
        storage_level_end = [inter(s.stop) for s in model_horizon_segments[1:end-1]]
        push!(storage_level_end, 0.5)
        storages[es].storage_level_start = storage_level_start
        storages[es].storage_level_end = storage_level_end
    end
end



