# -------------------------------------------------------------------------------------------------
# POMATO - Power Market Tool (C) 2018
# Current Version: Pre-Release-2018
# Created by Robert Mieth and Richard Weinhold
# Licensed under LGPL v3
#
# Language: Julia, v1.1.0 (required)
# ----------------------------------
#
# This file:
# Julia mutable struct definitions for efficient data handling
# -------------------------------------------------------------------------------------------------



mutable struct Line
    # Attributes
    index::Int
    name::String
    b::Float64
    capacity::Float64
    zone_i::Int
    zone_j::Int
    incidence::Vector{Float64}
    function Line(index::Int,
                  name::String,
                  b::Float64,
                  capacity::Float64,
                  zone_i::Int,
                  zone_j::Int,
                  incidence::Vector{Int})
        l = new()
        l.index = index
        l.name = name
        l.b = b
        l.capacity = capacity
        l.zone_i = zone_i
        l.zone_j = zone_j
        l.incidence = incidence
        return l
    end
end

mutable struct Contingency
    # Attributes
    name::String
    # Angle Formulation
    lines::Vector{Int}
    outages::Vector{Int}
    # PTDF Formulation
    ptdf::Array{Float64}
    ram::Vector{Float64}
    # Optional Attributes
    timestep::String
    function Contingency(name::String,
                         lines::Vector{Int},
                         outages::Vector{Int},
                         ptdf::Array{Float64, 2},
                         ram::Vector{Float64})
        c = new()
        c.name = name
        c.lines = lines
        c.outages = outages
        c.ptdf = ptdf
        c.ram = ram
        return c
    end
end

mutable struct Zone
    # Attributes
    index::Int
    name::Any
    demand::Array
    nodes::Array
    plants::Array
    res_plants::Array
    # Optional Attributes
    net_position::Array
    net_export::Array
    ntc::Array
    function Zone(index::Int,
                  name::Any,
                  demand::Array,
                  nodes::Array,
                  plants::Array,
                  res_plants::Array)
        z = new()
        z.index = index
        z.name = name
        z.demand = demand
        z.nodes = nodes
        z.plants = plants
        z.res_plants = res_plants
        return z
    end
end

mutable struct Heatarea
    # Attributes
    index::Int
    name::Any
    demand::Array
    plants::Array
    res_plants::Array
    function Heatarea(index::Int,
                      name::Any,
                      demand::Array,
                      plants::Array,
                      res_plants::Array)
        ha = new()
        ha.index = index
        ha.name = name
        ha.demand = demand
        ha.plants = plants
        ha.res_plants = res_plants
        return ha
    end
end

mutable struct Node
    # Attributes
    index::Int
    name::Any
    zone::Int
    slack::Bool
    demand::Array
    demand_rt::Array
    real_time::Bool
    plants::Array
    res_plants::Array
    # Optional Attributes
    net_export::Array
    slack_zone::Array
    demand_da::Array
    function Node(index::Int,
                  name::Any,
                  zone::Int,
                  demand_rt::Array,
                  slack::Bool,
                  plants::Array,
                  res_plants::Array)
        n = new()
        n.index = index
        n.name = name
        n.zone = zone
        n.demand_rt = demand_rt
        n.demand = demand_rt
        n.real_time = true

        n.slack = slack
        n.plants = plants
        n.res_plants = res_plants
        return n
    end
end

mutable struct Renewables
        index::Int
        name::Any
        g_max::Float64
        h_max::Float64
        mc_el::Float64
        mc_heat::Float64

        real_time::Bool
        mu::Array
        mu_rt::Array
        mu_da::Array

        mu_heat::Array
        mu_heat_rt::Array
        mu_heat_da::Array

        sigma_factor::Float64
        sigma::Array
        sigma_rt::Array
        sigma_da::Array

        sigma_heat::Array
        sigma_heat_rt::Array
        sigma_heat_da::Array
        node::Int
        plant_type::Any
        function Renewables(index::Int, name::Any,
                            g_max::Float64, h_max::Float64,
                            mc_el::Float64, mc_heat::Float64,
                            availability_rt::Array, node::Int, plant_type::Any)
            res = new()
            res.index = index
            res.g_max = g_max
            res.h_max = h_max
            res.mc_el = mc_el
            res.mc_heat = mc_heat
            res.name = name
            res.node = node
            res.plant_type = plant_type

            res.sigma_factor =  0.2
            
            res.mu_rt = availability_rt * g_max
            res.sigma_rt = res.mu_rt * res.sigma_factor

            res.mu_heat_rt = availability_rt * h_max
            res.sigma_heat_rt = res.mu_heat_rt * res.sigma_factor
            
            res.mu = res.mu_rt
            res.sigma = res.sigma_rt
            res.mu_heat = res.mu_heat_rt
            res.sigma_heat = res.sigma_heat_rt

            res.real_time = true
            return res
        end
    end


mutable struct Plant
    # Attributes
    index::Int
    name::Any
    node::Int
    mc_el::Float64
    mc_heat::Float64
    g_max::Float64
    d_max::Float64
    h_max::Float64
    eta::Float64
    plant_type::Any
    # Optional Attributes
    storage_capacity::Float64
    inflow::Array
    storage_start::Float64
    storage_end::Float64
    storage_level::Array
    
    function Plant(index::Int,
                   name::Any,
                   node::Int,
                   mc_el::Float64,
                   mc_heat::Float64,
                   eta::Float64,
                   g_max::Float64,
                   h_max::Float64,
                   plant_type::Any)
        p = new()
        p.index = index
        p.name = name
        p.node = node
        p.mc_el = mc_el
        p.mc_heat = mc_heat
        p.eta = eta
        p.g_max = g_max
        p.h_max = h_max
        p.plant_type = plant_type
        return p
    end
end

mutable struct DC_Line
    # Attributes
    index::Int
    name::Any
    node_i::Int
    node_j::Int
    capacity::Float64
    function DC_Line(index::Int,
                     name::Any,
                     node_i::Int,
                     node_j::Int,
                     capacity::Float64)
        l = new()
        l.index = index
        l.name = name
        l.node_i = node_i
        l.node_j = node_j
        l.capacity = capacity
        return l
    end
end

mutable struct Timestep
    index::Int
    name::String
    function Timestep(index::Int, name::String)
        t = new()
        t.index = index
        t.name = name
        return t
    end
end

mutable struct Data
    # Attributes
    nodes::Vector{Node}
    zones::Vector{Zone}
    heatareas::Vector{Heatarea}
    plants::Vector{Plant}
    renewables::Vector{Renewables}
    lines::Vector{Line}
    contingencies::Vector{Contingency}
    redispatch_contingencies::Vector{Contingency}
    dc_lines::Vector{DC_Line}
    t::Vector{Timestep}
    folders::Dict{String, String}

    function Data(nodes::Vector{Node}, zones::Vector{Zone},
                  heatareas::Vector{Heatarea}, plants::Vector{Plant},
                  renewables::Vector{Renewables}, lines::Vector{Line},
                  contingencies::Vector{Contingency}, redispatch_contingencies::Vector{Contingency},
                  dc_lines::Vector{DC_Line}, t::Vector{Timestep})
          d = new()
          d.nodes = nodes
          d.zones = zones
          d.heatareas = heatareas
          d.plants = plants
          d.renewables = renewables
          d.lines = lines
          d.contingencies = contingencies
          d.redispatch_contingencies = redispatch_contingencies
          d.dc_lines = dc_lines
          d.t = t
          return d
      end
end
