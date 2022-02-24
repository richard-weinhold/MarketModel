"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Definition of the POMATO and Result struct that house and facilitate the model process.
"""

mutable struct Result
	G::DataFrame
	H::DataFrame
	INJ::DataFrame
	F_DC::DataFrame
	EX::DataFrame
	D_es::DataFrame
	L_es::DataFrame
	D_hs::DataFrame
	L_hs::DataFrame
	D_ph::DataFrame
	INFEASIBILITY_H_POS::DataFrame
	INFEASIBILITY_H_NEG::DataFrame
	INFEASIBILITY_EL_POS::DataFrame
	INFEASIBILITY_EL_NEG::DataFrame
	INFEASIBILITY_ES::DataFrame
	EB_nodal::DataFrame
	EB_zonal::DataFrame
	CURT::DataFrame
	Alpha::DataFrame
	CC_LINE_MARGIN::DataFrame
	INFEASIBILITY_CC_LINES::DataFrame
	G_RES::DataFrame
	H_RES::DataFrame
	COST_G::DataFrame
	COST_H::DataFrame
	COST_EX::DataFrame
	COST_CURT::DataFrame
	COST_REDISPATCH::DataFrame
	COST_INFEASIBILITY_EL::DataFrame
	COST_INFEASIBILITY_H::DataFrame
	COST_INFEASIBILITY_ES::DataFrame
	misc_results::Dict
	options::Dict
	function Result()
		return new()
	end
end

mutable struct POMATO
    ### Main Attributes
    model::Model
    data::Data
    options::Dict
	result::Result

    ### Mappings and Sets
    n::NamedTuple{(
		:t, :zones, :nodes, :heatareas,
        :plants, :res, :dc, :lines, :contingencies,
        :he, :chp, :es, :hs, :ph, :alpha, :cc_res), Tuple{Vararg{Int, 16}}}

    ## Plant Mappings
    mapping::NamedTuple{(
		:slack, # slacks to 1:n_nodes
        :he, # 1:N.he index to 1:N.plants
        :chp, # 1:N.chp to 1:N.he
        :es, # 1:N.es to 1:N.plants
        :hs, # 1:N.hs to 1:N.he
        :ph, # 1:N.ph to 1:N.he
        :alpha, # 1:N.alpha to 1:N.he
        :cc_res, # map 1:cc_res to 1:n_res
		), Tuple{Vararg{Vector{Int}, 8}}}
		function POMATO()
			return new()
		end
end

function POMATO(model::Model, data::Data)

	m = POMATO()
	m.model = model
	m.data = data
	m.options = data.options

	## Plant Mappings
	# mapping heat index to G index
	mapping_he = findall(plant -> plant.h_max > 0, data.plants)

	m.mapping = (
		slack = findall(node -> node.slack, data.nodes),
		he = mapping_he,
		chp = findall(plant -> ((plant.h_max > 0)&(plant.g_max > 0)), data.plants[mapping_he]),
		es = findall(plant -> plant.plant_type in m.options["plant_types"]["es"], data.plants),
		hs = findall(plant -> plant.plant_type in m.options["plant_types"]["hs"], data.plants[mapping_he]),
		ph = findall(plant -> plant.plant_type in m.options["plant_types"]["ph"], data.plants[mapping_he]),
		alpha = findall(plant -> ((plant.g_max > m.options["chance_constrained"]["alpha_plants_mw"])&(plant.mc_el <= m.options["chance_constrained"]["alpha_plants_mc"])), data.plants),
		cc_res = findall(res_plants -> res_plants.g_max > m.options["chance_constrained"]["cc_res_mw"], data.renewables),
	)

	m.n = (t = size(data.t, 1),
		   zones = size(data.zones, 1),
		   nodes = size(data.nodes, 1),
		   heatareas = size(data.heatareas, 1),
		   plants = size(data.plants, 1),
		   res = size(data.renewables, 1),
		   dc = size(data.dc_lines, 1),
		   lines = size(data.lines, 1),
		   contingencies = size(data.contingencies, 1),
		   he = size(m.mapping.he, 1),
		   chp = size(m.mapping.chp, 1),
		   es = size(m.mapping.es, 1),
		   hs = size(m.mapping.hs, 1),
		   ph = size(m.mapping.ph, 1),
		   alpha = size(m.mapping.alpha, 1),
		   cc_res = size(m.mapping.cc_res, 1))
	return m
end

