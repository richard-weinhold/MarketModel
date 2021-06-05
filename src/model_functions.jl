

function add_variables_expressions!(pomato::POMATO)
	model = pomato.model
	n = pomato.n
	data = pomato.data
	options = pomato.options
	mapping = pomato.mapping

	@variable(model, G[1:n.t, 1:n.plants] >= 0) # El. power generation per plant p
	@variable(model, H[1:n.t, 1:n.he] >= 0) # Heat generation per plant p
	@variable(model, D_es[1:n.t, 1:n.es] >= 0) # El. demand of storage plants
	@variable(model, D_hs[1:n.t, 1:n.hs] >= 0) # El. demand of heat storage
	@variable(model, D_ph[1:n.t, 1:n.ph] >= 0) # El. demand of power to heat
	@variable(model, L_es[1:n.t, 1:n.es] >= 0) # Level of electricity storage
	@variable(model, L_hs[1:n.t, 1:n.hs] >= 0) # Level of heat storage

	@variable(model, EX[1:n.t, 1:n.zones, 1:n.zones] >= 0) # Commercial Exchanges between zones (row from, col to)
	@variable(model, INJ[1:n.t, 1:n.nodes]) # Net Injection at Node n
	@variable(model, F_DC_POS[1:n.t, 1:n.dc] >= 0) # Positive component of  dc line flow
	@variable(model, F_DC_NEG[1:n.t, 1:n.dc] >= 0) # Negative component of  dc line flow
	@expression(model, F_DC[t=1:n.t, dc=1:n.dc], F_DC_POS[t, dc] - F_DC_NEG[t, dc]) 

	if options["infeasibility"]["electricity"]["include"]
		infeasibility_bound = options["infeasibility"]["electricity"]["bound"]
	    @variable(model, 
			0 <= INFEASIBILITY_EL_NEG[1:n.t, 1:n.nodes] <= infeasibility_bound
		);
	    @variable(model, 
			0 <= INFEASIBILITY_EL_POS[1:n.t, 1:n.nodes] <= infeasibility_bound
		);
	else
	    @variable(model, INFEASIBILITY_EL_NEG[1:n.t, 1:n.nodes] == 0)
	    @variable(model, INFEASIBILITY_EL_POS[1:n.t, 1:n.nodes] == 0)
	end

	if options["infeasibility"]["heat"]["include"]
	    # Relaxing at high costs to avoid infeasibility in heat EB
		infeasibility_bound_heat = options["infeasibility"]["heat"]["bound"]
	    @variable(model, 
			0 <= INFEASIBILITY_H_NEG[1:n.t, 1:n.heatareas] <= infeasibility_bound_heat
		);
	    @variable(model, 
			0 <= INFEASIBILITY_H_POS[1:n.t, 1:n.heatareas] <= infeasibility_bound_heat
		);
	else
	    @expression(model, INFEASIBILITY_H_NEG[1:n.t, 1:n.heatareas], 0)
	    @expression(model, INFEASIBILITY_H_POS[1:n.t, 1:n.heatareas], 0)
	end

	@expression(model, INFEASIBILITY_EL_ZONAL_NEG[t=1:n.t, zone=1:n.zones],
	    sum(INFEASIBILITY_EL_NEG[t, node] for node in data.zones[zone].nodes));

	@expression(model, INFEASIBILITY_EL_ZONAL_POS[t=1:n.t, zone=1:n.zones],
	    sum(INFEASIBILITY_EL_POS[t, node] for node in data.zones[zone].nodes));

	# Expressions to create Nodal/Zonal Generation and Res Feedin
	@expression(model, G_Node[t=1:n.t, node=1:n.nodes],
		size(data.nodes[node].plants, 1) > 0 ? sum(G[t, plant] for plant in data.nodes[node].plants) : 0);

	@expression(model, D_Node[t=1:n.t, node=1:n.nodes],
		size(intersect(data.nodes[node].plants, mapping.ph), 1) > 0
	    ? sum(D_ph[t, findfirst(ph -> ph==plant, mapping.ph)] for plant in intersect(data.nodes[node].plants, mapping.ph))
	    : 0
		+ size(intersect(data.nodes[node].plants, mapping.es), 1) > 0
	    ? sum(D_es[t, findfirst(es -> es==plant, mapping.es)] for plant in intersect(data.nodes[node].plants, mapping.es))
	    : 0 );

	@expression(model, G_Zone[t=1:n.t, z=1:n.zones],
		sum(G_Node[t, node] for node in data.zones[z].nodes));

	@expression(model, D_Zone[t=1:n.t, z=1:n.zones],
		sum(D_Node[t, node] for node in data.zones[z].nodes));

	@expression(model, H_Heatarea[t=1:n.t, ha=1:n.heatareas],
		size(data.heatareas[ha].plants, 1) > 0
	    ? sum(H[t, findfirst(he -> he==plant, mapping.he)] for plant in data.heatareas[ha].plants)
	    : 0);

	@expression(model, D_Heatarea[t=1:n.t, ha=1:n.heatareas],
		size(intersect(data.heatareas[ha].plants, mapping.hs), 1) > 0
	    ? sum(D_hs[t, findfirst(hs -> hs==plant, mapping.hs)] for plant in intersect(data.heatareas[ha].plants, mapping.hs))
	    : 0);

	@expression(model, G_RES[t=1:n.t, res=1:n.res],
		GenericAffExpr{Float64, VariableRef}(data.renewables[res].mu[t]));
	@expression(model, H_RES[t=1:n.t, res=1:n.res],
		GenericAffExpr{Float64, VariableRef}(data.renewables[res].mu_heat[t]));

	@expression(model, RES_Node[t=1:n.t, node=1:n.nodes],
		size(data.nodes[node].res_plants, 1) > 0
	    ? GenericAffExpr{Float64, VariableRef}(sum(data.renewables[res].mu[t]*1. for res in data.nodes[node].res_plants))
	    : GenericAffExpr{Float64, VariableRef}(0));

	@expression(model, RES_Zone[t=1:n.t, z=1:n.zones],
	    sum(RES_Node[t, node] for node in data.zones[z].nodes));

	@expression(model, RES_Heatarea[t=1:n.t, ha=1:n.heatareas],
		size(data.heatareas[ha].res_plants, 1) > 0
	    ? sum(data.renewables[res].mu_heat[t] for res in data.heatareas[ha].res_plants)
	    : 0);

	add_cost_expressions!(pomato)
end

function add_cost_expressions!(pomato::POMATO)
	model = pomato.model
	n = pomato.n
	data = pomato.data
	options = pomato.options
	mapping = pomato.mapping

	G, H, EX = model[:G], model[:H], model[:EX]
	F_DC_POS, F_DC_NEG = model[:F_DC_POS], model[:F_DC_NEG]

	INFEASIBILITY_EL_POS, INFEASIBILITY_EL_NEG = model[:INFEASIBILITY_EL_POS], model[:INFEASIBILITY_EL_NEG]
	INFEASIBILITY_H_POS, INFEASIBILITY_H_NEG = model[:INFEASIBILITY_H_POS], model[:INFEASIBILITY_H_NEG]
	@expression(model, COST_G[t=1:n.t],
	    sum(G[t, p]*data.plants[p].mc_el for p in 1:n.plants));
	@expression(model, COST_H[t=1:n.t],
		n.he > 0 ? sum(H[t, he]*data.plants[mapping.he[he]].mc_heat for he in 1:n.he) : 0);
	if n.res > 0
		for t in 1:n.t
			add_to_expression!(COST_G[t], 
				sum(data.renewables[res].mu[t] * data.renewables[res].mc_el for res in 1:n.res)
			);
		end
	end
	if (n.res > 0)&(n.he > 0)
		for t in 1:n.t
			add_to_expression!(COST_H[t], 
				sum(data.renewables[res].mu_heat[t] * data.renewables[res].mc_heat for res in 1:n.res)
			);
		end
	end
	@expression(model, COST_EX[t=1:n.t], sum(EX[t, :, :]) + sum(F_DC_POS[t, :]) + sum(F_DC_NEG[t, :]));
	
	capacity_at_node(n) = 
    	((length(n.plants) > 0 ? sum(data.plants[p].g_max for p in n.plants) : 0 ))

	cap_threshold = quantile([capacity_at_node(n) for n in data.nodes], 0.90)
	dem_threshold = quantile([sum(n.demand) for n in data.nodes], 0.90)

	high_capacity_nodes = findall(n -> capacity_at_node(n) > cap_threshold, data.nodes)
	high_demand_nodes = findall(n -> sum(n.demand) > dem_threshold, data.nodes)

	# Prefer nodes with high (90% quantile) installed capacities for positive infeasibility
	# nodes with high demand for negative infeasibility
	@expression(model, COST_INFEASIBILITY_EL[t=1:n.t], 
		(sum(n in high_capacity_nodes ? 0.9*INFEASIBILITY_EL_POS[t, n] : 1.1*INFEASIBILITY_EL_POS[t, n] for n in 1:n.nodes)
		 + sum(n in high_demand_nodes ? 0.9*INFEASIBILITY_EL_NEG[t, n] : 1.1*INFEASIBILITY_EL_NEG[t, n] for n in 1:n.nodes)
		 )*options["infeasibility"]["electricity"]["cost"]);

	@expression(model, COST_INFEASIBILITY_H[t=1:n.t], (sum(INFEASIBILITY_H_POS[t, :])
		+ sum(INFEASIBILITY_H_NEG[t, :]))*options["infeasibility"]["heat"]["cost"]);
	@expression(model, COST_CURT[t=1:n.t], GenericAffExpr{Float64, VariableRef}(0));
	@expression(model, COST_REDISPATCH[t=1:n.t], GenericAffExpr{Float64, VariableRef}(0));
end

function add_objective!(pomato::POMATO)
	model = pomato.model
	# Objective Function based on Cost Expressions
	@objective(model, Min, sum(model[:COST_G]) + sum(model[:COST_H]) + sum(model[:COST_EX])
						   + sum(model[:COST_INFEASIBILITY_EL]) + sum(model[:COST_INFEASIBILITY_H])
						   + sum(model[:COST_CURT]) + sum(model[:COST_REDISPATCH]));
end

function add_result!(pomato::POMATO)
	pomato.result = Result(pomato)
end

function add_electricity_storage_constraints!(pomato::POMATO)
	model, n, mapping, data, options = pomato.model, pomato.n, pomato.mapping, pomato.data, pomato.options
	D_es, L_es, G = model[:D_es], model[:L_es], model[:G]

	@variable(model, Dump_Water[1:n.t, 1:n.es] >= 0);
	# Electricity Storage Equations
	storage_start(es) = data.plants[mapping.es[es]].storage_start*data.plants[mapping.es[es]].storage_capacity
	@constraint(model, [t=1:n.t, es=1:n.es],
			L_es[t, es]  == (t>1 ? L_es[t-1, es] : storage_start(es)) + data.plants[mapping.es[es]].inflow[t]  
			- G[t, mapping.es[es]] - Dump_Water[t, es] + data.plants[mapping.es[es]].eta*D_es[t, es])

	@constraint(model, [t=1:n.t],
		L_es[t, :] .<= [data.plants[mapping.es[es]].storage_capacity for es in 1:n.es])

	@constraint(model, [t=1:n.t],
		D_es[t, :] .<= [data.plants[mapping.es[es]].d_max for es in 1:n.es])

	# lower bound on storage level in last timestep
	storage_end(es) = data.plants[mapping.es[es]].storage_end*data.plants[mapping.es[es]].storage_capacity
	@constraint(model, L_es[n.t, :] .>= [storage_end(es) for es in 1:n.es]);
end

function add_electricity_generation_constraints!(pomato::POMATO)
	model, n, mapping, data, options = pomato.model, pomato.n, pomato.mapping, pomato.data, pomato.options
	# make Variable References Available
	G, INJ, EX, F_DC = model[:G], model[:INJ], model[:EX], model[:F_DC]
	INFEASIBILITY_EL_POS = model[:INFEASIBILITY_EL_POS]
	INFEASIBILITY_EL_NEG = model[:INFEASIBILITY_EL_NEG]
	# make Expression References Available
	INFEASIBILITY_EL_ZONAL_POS = model[:INFEASIBILITY_EL_ZONAL_POS]
	INFEASIBILITY_EL_ZONAL_NEG = model[:INFEASIBILITY_EL_ZONAL_NEG]
	G_Node, G_Zone = model[:G_Node], model[:G_Zone]
	D_Node, D_Zone = model[:D_Node], model[:D_Zone]
	RES_Node, RES_Zone = model[:RES_Node], model[:RES_Zone]

	# G Upper Bound
	@constraint(model, [t=1:n.t],
		G[t, :] .<= [data.plants[p].g_max for p in 1:n.plants])

	# DC Lines Constraints
	@constraint(model, [t=1:n.t],
	    F_DC[t, :] .<= [data.dc_lines[dc].capacity for dc in 1:n.dc])
	@constraint(model, [t=1:n.t],
	    -F_DC[t, :] .<= [data.dc_lines[dc].capacity for dc in 1:n.dc])

	# Balance Net Injections within Slacks Zones
	@constraint(model, [t=1:n.t, slack=mapping.slack],
	    0 == sum(INJ[t, n] for n in data.nodes[slack].slack_zone))
end

function add_electricity_energy_balance!(pomato::POMATO)
	model, n, data = pomato.model, pomato.n, pomato.data

	G_Zone, G_Node = model[:G_Zone], model[:G_Node]
	RES_Zone, RES_Node = model[:RES_Zone], model[:RES_Node]
	D_Zone, D_Node = model[:D_Zone], model[:D_Node]
	INFEASIBILITY_EL_POS = model[:INFEASIBILITY_EL_POS]
	INFEASIBILITY_EL_NEG = model[:INFEASIBILITY_EL_NEG]
	# make Expression References Available
	INFEASIBILITY_EL_ZONAL_POS = model[:INFEASIBILITY_EL_ZONAL_POS]
	INFEASIBILITY_EL_ZONAL_NEG = model[:INFEASIBILITY_EL_ZONAL_NEG]
	INJ = model[:INJ]
	F_DC = model[:F_DC]
	EX = model[:EX]

	# Create incidence matrix for dc-lines
	dc_incidence = spzeros(Int, n.nodes, n.dc)
	for dc in 1:n.dc
		dc_incidence[data.dc_lines[dc].node_i, dc] =  1
		dc_incidence[data.dc_lines[dc].node_j, dc] =  -1
	end

	# Zonal Energy Balance
	@constraint(model, EB_zonal[t=1:n.t, z=1:n.zones],
		data.zones[z].demand[t] ==
		+ data.zones[z].net_export[t]
		+ G_Zone[t, z] + RES_Zone[t, z] - D_Zone[t, z]
		- sum(EX[t, z, zz] for zz in 1:n.zones)
		+ sum(EX[t, zz, z] for zz in 1:n.zones)
		+ INFEASIBILITY_EL_ZONAL_POS[t, z] - INFEASIBILITY_EL_ZONAL_NEG[t, z]
		)

	# Nodal Energy Balance
	@constraint(model, EB_nodal[t=1:n.t, node=1:n.nodes],
		data.nodes[node].demand[t] ==
		+ data.nodes[node].net_export[t]
		+ G_Node[t, node] + RES_Node[t, node] - D_Node[t, node]
		- sum(dc_incidence[node, dc]*F_DC[t, dc] for dc in 1:n.dc)
		- INJ[t, node]
		+ INFEASIBILITY_EL_POS[t, node] - INFEASIBILITY_EL_NEG[t, node]
		)
end

function add_heat_generation_constraints!(pomato::POMATO)
	model, n, mapping, data, options = pomato.model, pomato.n, pomato.mapping, pomato.data, pomato.options
	chp_efficiency = ("chp_efficiency" in keys(options["parameters"]) ? options["parameters"]["chp_efficiency"] : 0.1)
	storage_start = options["parameters"]["storage_start"]

	# make Variable References Available
	G, H = model[:G], model[:H]
	D_hs, D_ph = model[:D_hs], model[:D_ph]
	L_hs = model[:L_hs]
	INFEASIBILITY_H_POS = model[:INFEASIBILITY_H_POS]
	INFEASIBILITY_H_NEG = model[:INFEASIBILITY_H_NEG]
	# make Expression References Available
	H_Heatarea = model[:H_Heatarea]
	D_Heatarea = model[:D_Heatarea]
	RES_Heatarea = model[:RES_Heatarea]

	# H upper bound
	@constraint(model, [t=1:n.t],
		H[t, :] .<= [data.plants[mapping.he[he]].h_max for he in 1:n.he])

	# CHP plants
	@constraint(model, [t=1:n.t, chp in 1:n.chp],
	    G[t, mapping.he[mapping.chp[chp]]] >= ((data.plants[mapping.he[mapping.chp[chp]]].g_max*(1 - chp_efficiency))
			/ data.plants[mapping.he[mapping.chp[chp]]].h_max) * H[t, mapping.chp[chp]])

	@constraint(model, [t=1:n.t, chp in 1:n.chp],
	    G[t, mapping.he[mapping.chp[chp]]] <= data.plants[mapping.he[mapping.chp[chp]]].g_max
									  * (1 - (chp_efficiency * H[t, mapping.chp[chp]] / data.plants[mapping.he[mapping.chp[chp]]].h_max)))
	# Power To Heat
	@constraint(model, [t=1:n.t, ph=1:n.ph],
		D_ph[t, ph] == H[t, mapping.ph[ph]] / data.plants[mapping.he[mapping.ph[ph]]].eta)
	@constraint(model, [t=1:n.t, ph=1:n.ph],
		G[t, ph] == 0)

	# Heat Storage Equations
	@constraint(model, [t=1:n.t, hs=1:n.hs],
	    L_hs[t, hs] ==  (t>1 ? data.plants[mapping.he[mapping.hs[hs]]].eta*L_hs[t-1, hs] : storage_start*data.plants[mapping.he[mapping.hs[hs]]].storage_capacity)
	                   - H[t, mapping.hs[hs]]
	                   + D_hs[t, hs])

	@constraint(model, [t=1:n.t],
	    L_hs[t, :] .<= [data.plants[mapping.he[mapping.hs[hs]]].storage_capacity for hs in 1:n.hs])
	@constraint(model, [t=1:n.t],
	    D_hs[t, :] .<= [data.plants[mapping.he[mapping.hs[hs]]].h_max for hs in 1:n.hs])

	# Heat Energy Balance
	@constraint(model, EB_Heat[t=1:n.t],
	    [data.heatareas[ha].demand[t] for ha in 1:n.heatareas] .==
	    .+ H_Heatarea[t, :] .+ RES_Heatarea[t, :] .- D_Heatarea[t, :]
	    .+ INFEASIBILITY_H_POS[t, :] .- INFEASIBILITY_H_NEG[t, :])
end

function add_curtailment_constraints!(pomato::POMATO, zones::Vector{String}, curt_market::Array{Float64, 2})
	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	G_RES = model[:G_RES]
	RES_Node = model[:RES_Node]
	RES_Zone = model[:RES_Zone]
	redispatch_zones_nodes = vcat([data.zones[z].nodes for z in findall(z -> z.name in zones, data.zones)]...)
	res_in_zone = findall(r -> r.node in redispatch_zones_nodes, data.renewables)

	@variable(model, CURT_REDISPATCH[1:n.t, res_in_zone] >= 0)
	@expression(model, CURT[t=1:n.t, res=1:n.res], GenericAffExpr{Float64, VariableRef}(0));
	for res in 1:n.res, t in 1:n.t
		add_to_expression!(CURT[t, res], res in res_in_zone ? CURT_REDISPATCH[t, res] : curt_market[t, res])
	end
	@constraint(model, MinCurt[t=1:n.t, res=res_in_zone],
		CURT[t, res] >= curt_market[t, res])
	@constraint(model, MaxCurt[t=1:n.t, res=res_in_zone],
		CURT[t, res] <= G_RES[t, res])
	for t in 1:n.t
		for res in 1:n.res
			add_to_expression!(G_RES[t, res],  -CURT[t, res])
		 	add_to_expression!(RES_Node[t, data.renewables[res].node], -CURT[t, res])
		 	add_to_expression!(RES_Zone[t, data.nodes[data.renewables[res].node].zone], -CURT[t, res])
	 	end
	 end
 	# @constraint(model, [t=1:n.t, res=res_in_zone], CURT[t, res] <= 0.2*data.renewables[res].mu[t])
	COST_CURT = model[:COST_CURT]
	for t in 1:n.t
		add_to_expression!(COST_CURT[t], sum(CURT[t, :])*pomato.options["curtailment"]["cost"])
	end
end

function add_curtailment_constraints!(pomato::POMATO)
	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	@variable(model, CURT[1:n.t, 1:n.res] >= 0)

	@constraint(model, MaxCurt[t=1:n.t, res=1:n.res],
		CURT[t, res] <= data.renewables[res].mu[t])

	G_RES = model[:G_RES]
	RES_Node = model[:RES_Node]
	RES_Zone = model[:RES_Zone]
	for t in 1:n.t
		for res in 1:n.res
			add_to_expression!(G_RES[t, res],  -CURT[t, res])
		 	add_to_expression!(RES_Node[t, data.renewables[res].node], -CURT[t, res])
		 	add_to_expression!(RES_Zone[t, data.nodes[data.renewables[res].node].zone], -CURT[t, res])
	 	end
	 end
	# @constraint(model, [t=1:n.t], CURT[t, :] .<= [0.2*res.mu[t] for res in data.renewables])
	COST_CURT = model[:COST_CURT]
	for t in 1:n.t
		add_to_expression!(COST_CURT[t], sum(CURT[t, :])*pomato.options["curtailment"]["cost"])
	end
end

function add_dclf_angle_constraints!(pomato::POMATO)
	add_dclf_angle_constraints!(pomato, collect(1:pomato.n.lines))
end

function add_dclf_angle_constraints!(pomato::POMATO, line_subset::Vector{Int})
	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	INJ = model[:INJ]
	
	@variable(model, THETA[1:n.t, 1:n.nodes, 1:n.contingencies]);
	@variable(model, F[1:n.t, 1:n.lines, 1:n.contingencies]);
	
	@constraint(model, [t=1:n.t, line=1:n.lines, contingency=1:n.contingencies], 
		F[t, line, contingency] == data.lines[line].b * sum((line in data.contingencies[contingency].outages ? 0 : data.lines[line].incidence[node]) * THETA[t, node, contingency] for node in 1:n.nodes));
	
	@constraint(model, [t=1:n.t, node=1:n.nodes, contingency=1:n.contingencies], 
		INJ[t, node] == sum(data.lines[line].incidence[node]*F[t, line, contingency] for line in 1:n.lines));

	@constraint(model, [t=1:n.t, slack=mapping.slack, contingency=1:n.contingencies], 0 == THETA[t, slack, contingency]);

	for contingency in 1:n.contingencies
		for line in intersect(data.contingencies[contingency].lines, line_subset)
			@constraint(model, [t=1:n.t], -data.lines[line].capacity <= F[t, line, contingency] <= data.lines[line].capacity);
		end
		for line in intersect(data.contingencies[contingency].outages, line_subset)
			@constraint(model, [t=1:n.t],  F[t, line, contingency] == 0);
		end
	end
end

function add_dclf_ptdf_constraints!(pomato::POMATO, line_subset::Vector{Int})
	model, n, data = pomato.model, pomato.n, pomato.data
	INJ = model[:INJ]
	tmp_ptdf = []
	tmp_capacity = []
	for contingency in data.contingencies
		line_indices = findall(l -> l in line_subset, contingency.lines)
		push!(tmp_ptdf, contingency.ptdf[line_indices, :])
		push!(tmp_capacity, contingency.ram[line_indices])
	end
	ptdf = vcat(tmp_ptdf...)
	capacity = vcat(tmp_capacity...)
	
	@variable(model, F_pos[1:n.t, 1:size(ptdf, 1)] >= 0)
	@variable(model, F_neg[1:n.t, 1:size(ptdf, 1)] >= 0)
	@constraint(model, [t=1:n.t], ptdf * INJ[t, :] .== F_pos[t, :] .- F_neg[t, :]);
	@constraint(model, [t=1:n.t], F_pos[t, :] .<= capacity);
	@constraint(model, [t=1:n.t], F_neg[t, :] .<= capacity);
end

function add_dclf_ptdf_constraints!(pomato::POMATO)
	model, n, data = pomato.model, pomato.n, pomato.data
	INJ = model[:INJ]
	ptdf = vcat([contingency.ptdf for contingency in data.contingencies]...)
	capacity = vcat([contingency.ram for contingency in data.contingencies]...)
	
	@variable(model, F_pos[1:n.t, 1:size(ptdf, 1)] >= 0)
	@variable(model, F_neg[1:n.t, 1:size(ptdf, 1)] >= 0)
	@constraint(model, [t=1:n.t], ptdf * INJ[t, :] .== F_pos[t, :] .- F_neg[t, :]);
	@constraint(model, [t=1:n.t], F_pos[t, :] .<= capacity);
	@constraint(model, [t=1:n.t], F_neg[t, :] .<= capacity);
end

function add_flowbased_constraints!(pomato::POMATO)
	model, n, data = pomato.model, pomato.n, pomato.data
	EX = model[:EX]
	# zonal PTDF * NEX <= Capacity
	# For FB Domain PTDF contains upper and lower bounds and it is time dependant
	# For the Zonal PTDF Bounds are symmetrical, therefore we need two constraints
	contingency_t(t) = filter(c -> c.timestep == data.t[t].name, data.contingencies)
	if any(isdefined(contingency, :timestep) for contingency in data.contingencies)
		@info("Adding FB Domain for each timestep")
		@constraint(model, [t=1:n.t], 
			vcat([contingency.ptdf for contingency in contingency_t(t)]...) 
			* sum(EX[t, :, zz] - EX[t, zz, :] for zz in 1:n.zones) 
			.<= vcat([contingency.ram for contingency in contingency_t(t)]...));
	else
		@info("Adding zonal PTDF")
		zonal_ptdf = vcat([contingency.ptdf for contingency in data.contingencies]...)
		ram = vcat([contingency.ram for contingency in data.contingencies]...)
		for t in 1:n.t
			nex = sum(EX[t, :, zz] - EX[t, zz, :] for zz in 1:n.zones)
			@constraint(model, zonal_ptdf * nex .<= ram);
			@constraint(model, -zonal_ptdf * nex .<= ram);
		end
	end
end

function add_chance_constrained_flowbased_constraints!(pomato::POMATO)

	model, n, mapping, data, options = pomato.model, pomato.n, pomato.mapping, pomato.data, pomato.options;
	# mapping CC Res to Nodes
	cc_res_to_zone = spzeros(Int8, n.zones, n.res)
	for res in data.renewables[mapping.cc_res]
		cc_res_to_zone[data.nodes[data.renewables[mapping.cc_res][1].node].zone, res.index] = Int8(1)
	end
	cc_res_to_zone = cc_res_to_zone[:, mapping.cc_res]

	# Distribution Paramters
	epsilon = 0.05
	z = quantile(Normal(0,1), 1-epsilon);
	sigma = hcat([res.sigma for res in data.renewables[mapping.cc_res]]...);
	# sigma = zeros(n.t, n_cc_res)
	Epsilon = [sparse(diagm(0 => (sigma[t, :].^2))) for t in 1:n.t];
	Epsilon_root = [Epsilon[t]^(1/2) for t in 1:n.t];
	S = [sqrt(sum(Epsilon[t])) for t in 1:n.t];

	alpha_plants(node) = findall(x -> x in intersect(mapping.alpha, data.nodes[node].plants), mapping.alpha);
	if options["chance_constrained"]["fixed_alpha"]
		@expression(model, Alpha[t=1:n.t, alpha=1:n.alpha],
				data.plants[mapping.alpha[alpha]].g_max/(sum(pp.g_max for pp in data.plants[mapping.alpha])));
		@expression(model, Alpha_Nodes[t=1:n.t, node=1:n.nodes],
			size(alpha_plants(node), 1) > 0 ?
			sum(Alpha[t, alpha] for alpha in alpha_plants(node)) : 0);
	else
		@variable(model, Alpha[1:n.t, 1:n.alpha] >= 0);
		@expression(model, Alpha_Nodes[t=1:n.t, node=1:n.nodes],
			size(alpha_plants(node), 1) > 0 ?
			sum(Alpha[t, alpha] for alpha in alpha_plants(node)) : 0);
		Alpha_Nodes = convert(Array{GenericAffExpr{Float64, VariableRef},2}, Alpha_Nodes)
	end
				
	@expression(model, Alpha_Zones[t=1:n.t, z=1:n.zones],
		sum(Alpha_Nodes[t, node] for node in data.zones[z].nodes));
	G, INJ = model[:G], model[:INJ];
	@constraint(model, [t=1:n.t], sum(Alpha[t, :]) == 1);
	@constraint(model, [t=1:n.t, alpha=1:n.alpha], G[t, mapping.alpha[alpha]] + z*Alpha[t,alpha]*S[t] <= data.plants[mapping.alpha[alpha]].g_max);
	@constraint(model, [t=1:n.t, alpha=1:n.alpha], -G[t, mapping.alpha[alpha]] + z*Alpha[t,alpha]*S[t] <= 0);

	@info("Adding load flow constaints... ")
	contingency_t(t) = filter(c -> c.timestep == data.t[t].name, data.contingencies)
	contingency_idx_t(t) = findall(c -> c.timestep == data.t[t].name, data.contingencies)
	ptdf_t(t) = vcat([contingency.ptdf for contingency in contingency_t(t)]...) 
	ram_t(t) = vcat([contingency.ram for contingency in contingency_t(t)]...) 
	@variable(model, 
		CC_LINE_MARGIN[t=1:n.t, co=contingency_idx_t(t), cb=1:length(data.contingencies[co].lines)] >= 0
	);
	CC_LINE_MARGIN_t(t) = [CC_LINE_MARGIN[t, co, cb] for co in contingency_idx_t(t) for cb in 1:length(data.contingencies[co].lines)]

	EX = model[:EX]
	for t in 1:n.t
		ptdf = ptdf_t(t)
		ram = ram_t(t)
		B = cc_res_to_zone - hcat([Alpha_Zones[t,:] for i in 1:size(cc_res_to_zone, 2)]...)
		alpha_loadflow = create_alpha_loadflow_constraint(ptdf, B, cc_res_to_zone)
		@constraint(model, [cb=1:size(ptdf, 1)], 
			vec(vcat(CC_LINE_MARGIN_t(t)[cb], (alpha_loadflow[cb, :]'*Epsilon_root[t])')) in SecondOrderCone()
		);
		@constraint(model, 
			Array(ptdf * sum(EX[t, :, zz] - EX[t, zz, :] for zz in 1:n.zones)/z .+ CC_LINE_MARGIN_t(t)) .<= ram/z
		);
	end
end

function add_ntc_constraints!(pomato::POMATO)
	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	EX = model[:EX]
	@constraint(model, [t=1:n.t, z=1:n.zones, zz=1:n.zones], 
		EX[t, z, zz] <= data.zones[z].ntc[zz]
	);
end

function add_ntc_constraints!(pomato::POMATO, zones::Vector{Int})

	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	EX = model[:EX]
	@info("Including NTC Constraints for: "*join([data.zones[z].name*", " for z in zones])[1:end-2])
	for non_fb_zone in zones
		@constraint(model, [t=1:n.t, zz=1:n.zones], 
			EX[t, non_fb_zone, zz] <= data.zones[non_fb_zone].ntc[zz]
		);
		@constraint(model, [t=1:n.t, z=1:n.zones], 
			EX[t, z, non_fb_zone] <= data.zones[z].ntc[non_fb_zone]
		);
	end
end

function add_net_position_constraints!(pomato::POMATO)
	model, n, mapping, data = pomato.model, pomato.n, pomato.mapping, pomato.data
	EX = model[:EX]
	nex_zones = [z.index for z in data.zones if any(z.net_position .!== missing)]
	@info("Including NEX Constraints for: "*join([data.zones[z].name*", " for z in nex_zones])[1:end-2])
	# # Applies to nodal model for basecase calculation:
	@constraint(model, [t=1:n.t, z=nex_zones], 
		sum(EX[t, z, zz] - EX[t, zz, z] for zz in 1:n.zones)
		<= data.zones[z].net_position[t] + 0.1*abs(data.zones[z].net_position[t])
	);
	@constraint(model, [t=1:n.t, z=nex_zones], 
		sum(EX[t, z, zz] - EX[t, zz, z] for zz in 1:n.zones)
		>= data.zones[z].net_position[t] - 0.1*abs(data.zones[z].net_position[t])
	);
end

function create_alpha_loadflow_constraint(
	ptdf::Array{Float64, 2},
	B::Union{Array{GenericAffExpr{Float64, VariableRef},2}, Array{Real,2}, Array{Float64,2}},
	cc_res_incidence::SparseMatrixCSC{Int8,Int64})

	nonzero_indices = [findall(x -> x != GenericAffExpr{Float64, VariableRef}(0), B[:, i]) for i in 1:size(cc_res_incidence, 2)]
	alpha_loadflow = zeros(GenericAffExpr{Float64, VariableRef}, 
						   size(ptdf, 1), size(cc_res_incidence, 2))
	# Create Segments
	n_hat = size(alpha_loadflow, 2)
	S = ceil(n_hat/Threads.nthreads())
	ranges = Array{UnitRange{Int},1}()
	for i in 1:Threads.nthreads()
		push!(ranges, range(Int((i-1)*S + 1), stop=Int(min(n_hat, i*S))))
	end
	Threads.@threads for r in ranges
		for i in r
			@inbounds for cb in 1:size(ptdf, 1)
	 			alpha_loadflow[cb, i] = sum(ptdf[cb, j]*B[j,i] for j in nonzero_indices[i])
	   		end
		end
	end
	return alpha_loadflow
end

function add_chance_constraints!(pomato::POMATO)
	model, n, mapping, data, options = pomato.model, pomato.n, pomato.mapping, pomato.data, pomato.options
	@info("Creating Chance Constraints for $(n.cc_res) RES generators which are larger than "
		  *"$(options["chance_constrained"]["cc_res_mw"]) MW")
	
	if options["chance_constrained"]["fixed_alpha"]
		@info("Reserve allocation with fixed ratio for $(n.alpha) power plants which are larger than "
			  *"$(options["chance_constrained"]["alpha_plants_mw"]) MW")
	else
		@info("Reserve is freely allocated for $(n.alpha) power plants which are larger than "
			  *"$(options["chance_constrained"]["alpha_plants_mw"]) MW")
	end
	# mapping CC Res to Nodes
	cc_res_to_node = spzeros(Int8, n.nodes, n.res)
	for res in data.renewables[mapping.cc_res]
		cc_res_to_node[res.node, res.index] = Int8(1)
	end
	cc_res_to_node = cc_res_to_node[:, mapping.cc_res]

	# Distribution Paramters
	epsilon = 0.01
	z = quantile(Normal(0,1), 1-epsilon)
	sigma = hcat([res.sigma for res in data.renewables[mapping.cc_res]]...)
	# sigma = zeros(n.t, n_cc_res)
	Epsilon = [sparse(diagm(0 => (sigma[t, :].^2))) for t in 1:n.t]
	Epsilon_root = [Epsilon[t]^(1/2) for t in 1:n.t]
	S = [sqrt(sum(Epsilon[t])) for t in 1:n.t]

	plants_node_alpha(node) = findall(x -> x in intersect(mapping.alpha, data.nodes[node].plants), mapping.alpha)
	if options["chance_constrained"]["fixed_alpha"]
		@expression(model, Alpha[t=1:n.t, alpha=1:n.alpha],
				data.plants[mapping.alpha[alpha]].g_max/(sum(pp.g_max for pp in data.plants[mapping.alpha])))
		
		@expression(model, Alpha_Nodes[t=1:n.t, node=1:n.nodes],
			size(plants_node_alpha(node), 1) > 0 ? sum(Alpha[t, alpha] for alpha in plants_node_alpha(node)) : 0
		);
	else
		@variable(model, Alpha[1:n.t, 1:n.alpha] >= 0)
		@expression(model, Alpha_Nodes[t=1:n.t, node=1:n.nodes],
			size(plants_node_alpha(node), 1) > 0 ? sum(Alpha[t, alpha] for alpha in plants_node_alpha(node)) : 0
		);
		Alpha_Nodes = convert(Array{GenericAffExpr{Float64, VariableRef},2}, Alpha_Nodes)
	end

	G, INJ = model[:G], model[:INJ]
	@constraint(model, [t=1:n.t], sum(Alpha[t, :]) == 1);
	@constraint(model, [t=1:n.t, alpha=1:n.alpha], 
		G[t, mapping.alpha[alpha]] + z*Alpha[t,alpha]*S[t] <= data.plants[mapping.alpha[alpha]].g_max
	);
	@constraint(model, [t=1:n.t, alpha=1:n.alpha], 
		-G[t, mapping.alpha[alpha]] + z*Alpha[t,alpha]*S[t] <= 0
	);

	@info("Adding load flow constaints... ")
	ptdf = vcat([contingency.ptdf for contingency in data.contingencies]...)
	capacity = vcat([contingency.ram for contingency in data.contingencies]...)
	@variable(model, CC_LINE_MARGIN[t=1:n.t, co=1:n.contingencies, cb=1:length(data.contingencies[co].lines)] >= 0)
	CC_LINE_MARGIN_t(t) = [CC_LINE_MARGIN[t, co, cb] for co in 1:n.contingencies for cb in 1:length(data.contingencies[co].lines)]

	for t in 1:n.t
		@info("Load Alpha for t = $(t) using $(Threads.nthreads()) Threads")
		B = cc_res_to_node - hcat([Alpha_Nodes[t,:] for i in 1:size(cc_res_to_node, 2)]...)
		# Create Segments
		alpha_loadflow = create_alpha_loadflow_constraint(ptdf, B, cc_res_to_node)
		@constraint(model, [cb=1:size(ptdf, 1)], 
			vec(vcat(CC_LINE_MARGIN_t(t)[cb], (alpha_loadflow[cb, :]'*Epsilon_root[t])')) in SecondOrderCone()
		);
	end
	@info("PTDF constraints... ")
 	@variable(model, F_pos[1:n.t, 1:size(ptdf, 1)] >= 0)
	@variable(model, F_neg[1:n.t, 1:size(ptdf, 1)] >= 0)

	@constraint(model, [t=1:n.t], ptdf * INJ[t, :] .== F_pos[t, :] .- F_neg[t, :]);
	@constraint(model, [t=1:n.t], F_pos[t, :] .+ z*CC_LINE_MARGIN_t(t) .<= capacity);
	@constraint(model, [t=1:n.t], F_neg[t, :] .+ z*CC_LINE_MARGIN_t(t) .<= capacity);
end

function redispatch_model!(pomato::POMATO, market_model_results::Dict, redispatch_zones::Vector{String})
	#### Previous solution and data
	n = pomato.n
	data = pomato.data
	model = pomato.model
	mapping = pomato.mapping
	g_market = market_model_results["g_market"]
	d_es_market = market_model_results["d_es_market"]
	d_ph_market = market_model_results["d_ph_market"]
	infeas_neg_market, infeas_pos_market = market_model_results["infeas_neg_market"], market_model_results["infeas_pos_market"]

	redispatch_zones_nodes = vcat([data.zones[z].nodes for z in findall(z -> z.name in redispatch_zones, data.zones)]...)
	redispatch_zones_plants = setdiff(findall(p -> p.node in redispatch_zones_nodes, data.plants), union(pomato.mapping.es, mapping.ph))

	## Build Redispatch model
	@variable(model, G_redispatch[1:n.t, redispatch_zones_plants] >= 0) # El. power generation per plant p
	@variable(model, G_redispatch_pos[1:n.t, redispatch_zones_plants] >= 0) # El. power generation per plant p
	@variable(model, G_redispatch_neg[1:n.t, redispatch_zones_plants] >= 0) # El. power generation per plant p

	@variable(model, EX[1:n.t, 1:n.zones, 1:n.zones] >= 0) # Commercial Exchanges between zones (row from, col to)
	@variable(model, INJ[1:n.t, 1:n.nodes]) # Net Injection at Node n
	@variable(model, F_DC_POS[1:n.t, 1:n.dc] >= 0) # Positive component of  dc line flow
	@variable(model, F_DC_NEG[1:n.t, 1:n.dc] >= 0) # Negative component of  dc line flow
	@expression(model, F_DC[t=1:n.t, dc=1:n.dc], F_DC_POS[t, dc] - F_DC_NEG[t, dc]) 


	@variable(model, 
		0 <= INFEAS_NEG[1:n.t, 1:n.nodes] <= pomato.options["infeasibility"]["electricity"]["bound"]
	);
	@variable(model, 
		0 <= INFEAS_POS[1:n.t, 1:n.nodes] <= pomato.options["infeasibility"]["electricity"]["bound"]
	);
	#
	@expression(model, 
		INFEASIBILITY_EL_NEG[t=1:n.t, node=1:n.nodes], GenericAffExpr{Float64, VariableRef}(0)
	);
	@expression(model, 
		INFEASIBILITY_EL_POS[t=1:n.t, node=1:n.nodes],  GenericAffExpr{Float64, VariableRef}(0)
	);
	@expression(model, 
		H[1:n.t, 1:n.he],  GenericAffExpr{Float64, VariableRef}(0)
	);
	@expression(model, 
		INFEASIBILITY_H_POS[1:n.t, 1:n.heatareas], GenericAffExpr{Float64, VariableRef}(0)
	);
	@expression(model, 
		INFEASIBILITY_H_NEG[1:n.t, 1:n.heatareas], GenericAffExpr{Float64, VariableRef}(0)
	);

	for node in 1:n.nodes, t in 1:n.t
		add_to_expression!(INFEASIBILITY_EL_NEG[t, node], INFEAS_NEG[t, node] + infeas_neg_market[t, node])
		add_to_expression!(INFEASIBILITY_EL_POS[t, node], INFEAS_POS[t, node] + infeas_pos_market[t, node])
	end

	@expression(model, 
		INFEASIBILITY_EL_ZONAL_NEG[t=1:n.t, z=1:n.zones],  GenericAffExpr{Float64, VariableRef}(0)
	);
	@expression(model, 
		INFEASIBILITY_EL_ZONAL_POS[t=1:n.t, z=1:n.zones],  GenericAffExpr{Float64, VariableRef}(0)
	);

	for node in 1:n.nodes, t in 1:n.t
		add_to_expression!(INFEASIBILITY_EL_ZONAL_NEG[t, data.nodes[node].zone], 
						   INFEASIBILITY_EL_NEG[t, node])
		add_to_expression!(INFEASIBILITY_EL_ZONAL_POS[t, data.nodes[node].zone], 
						   INFEASIBILITY_EL_POS[t, node])
	end

	# Expressions to create Nodal/Zonal Generation and Res Feedin
	@expression(model, G[t=1:n.t, p=1:n.plants], GenericAffExpr{Float64, VariableRef}(0));
	for p in 1:n.plants, t in 1:n.t
		add_to_expression!(G[t, p], p in redispatch_zones_plants ? G_redispatch[t,p] : g_market[t, p])
	end

	@expression(model, G_Node[t=1:n.t, node=1:n.nodes],
		size(data.nodes[node].plants, 1) > 0 ? sum(G[t, plant] for plant in data.nodes[node].plants) : 0);

	@expression(model, D_es[t=1:n.t, es=1:n.es], d_es_market[t, es])
	@expression(model, D_ph[t=1:n.t, ph=1:n.ph], d_ph_market[t, ph])

	@expression(model, D_Node[t=1:n.t, node=1:n.nodes],
		size(intersect(data.nodes[node].plants, mapping.ph), 1) > 0
	    ? sum(D_ph[t, findfirst(ph -> ph==plant, mapping.ph)] for plant in intersect(data.nodes[node].plants, mapping.ph))
	    : 0
		+ size(intersect(data.nodes[node].plants, mapping.es), 1) > 0
	    ? sum(D_es[t, findfirst(es -> es==plant, mapping.es)] for plant in intersect(data.nodes[node].plants, mapping.es))
	    : 0 );

	@expression(model, G_Zone[t=1:n.t, z=1:n.zones],
		sum(G_Node[t, node] for node in data.zones[z].nodes)
	);

	@expression(model, D_Zone[t=1:n.t, z=1:n.zones],
		sum(D_Node[t, node] for node in data.zones[z].nodes)
	);

	@expression(model, G_RES[t=1:n.t, res=1:n.res], 
		GenericAffExpr{Float64, VariableRef}(data.renewables[res].mu[t])
	);

	@expression(model, RES_Node[t=1:n.t, node=1:n.nodes],
		size(data.nodes[node].res_plants, 1) > 0
		? GenericAffExpr{Float64, VariableRef}(sum(data.renewables[res].mu[t] for res in data.nodes[node].res_plants))
		: GenericAffExpr{Float64, VariableRef}(0));

	@expression(model, RES_Zone[t=1:n.t, z=1:n.zones],
		sum(RES_Node[t, node] for node in data.zones[z].nodes));

	add_cost_expressions!(pomato);
	for t in 1:n.t
		add_to_expression!(model[:COST_REDISPATCH][t], 
			(sum(G_redispatch_pos[t, :]*pomato.options["redispatch"]["cost"]) + sum(G_redispatch_neg[t, :]*pomato.options["redispatch"]["cost"]))
		);
	end
	# G Upper Bound
	@constraint(model, [t=1:n.t, p=redispatch_zones_plants],
		G_redispatch[t, p] <= data.plants[p].g_max)
	# @constraint(model, [t=1:n.t, p=redispatch_zones_plants],
	# 	G_redispatch_neg[t, p] <= data.plants[p].g_max)
	# @constraint(model, [t=1:n.t, p=redispatch_zones_plants],
	# 	G_redispatch_pos[t, p] <= data.plants[p].g_max)
	@constraint(model, [t=1:n.t, p=redispatch_zones_plants],
		G_redispatch[t, p] - g_market[t, p] == G_redispatch_pos[t, p] - G_redispatch_neg[t, p])
	
	@constraint(model, [t=1:n.t],
		F_DC[t, :] .<= [data.dc_lines[dc].capacity for dc in 1:n.dc])
	@constraint(model, [t=1:n.t],
		-F_DC[t, :] .<= [data.dc_lines[dc].capacity for dc in 1:n.dc])

	# Balance Net Injections within Slacks Zones
	@constraint(model, [t=1:n.t, slack=mapping.slack],
		0 == sum(INJ[t, n] for n in data.nodes[slack].slack_zone))
		
	redispatch_zones_idx = map(name -> findfirst(zone -> zone.name == name, data.zones), redispatch_zones)  
	redispatch_line_subset = findall(line -> (line.zone_i in redispatch_zones_idx)&(line.zone_j in redispatch_zones_idx), data.lines)
	@info("$(length(redispatch_line_subset)) lines part of the redispatch network")
	add_dclf_angle_constraints!(pomato, redispatch_line_subset)
	add_objective!(pomato);
end
