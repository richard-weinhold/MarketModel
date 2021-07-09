"""
POMATO - Power Market Tool (C) 2021
Current Version: 0.4
Created by Richard Weinhold and Robert Mieth
Licensed under LGPL v3

Language: Julia, v1.5
----------------------------------

This file:
Gurobi specific functionality.
check_infeasibility allows to identify constraints that cause infeasibility when using Gurobi.

"""

function check_infeasibility(pomato::POMATO)
	global optimizer
	@info("Using Gurobi compute_conflict to find reason for infeasibility.")
	if string(optimizer) == "Gurobi.Optimizer"
		global optimizer_package
		model = pomato.model
		MOI.compute_conflict!(model.moi_backend.optimizer.model)
		for constraint_types in list_of_constraint_types(model)
			out = filter(x -> MOI.ConflictParticipationStatusCode(0) != 
								MOI.get(model.moi_backend, MOI.ConstraintConflictStatus(), x.index), 
							MarketModel.all_constraints(model, constraint_types[1], constraint_types[2]));
			
			constraint_identifier = join(constraint_types, "_")
			@warn("$(length(out)) constraints in conflic of type $constraint_identifier")
			CSV.write(pomato.data.folders["result_dir"]*"/"*constraint_identifier*".csv", 
						DataFrame([out[i] for i in axes(out)]))
		end
		@warn("Saved all constraints in conflict to $(pomato.data.folders["result_dir"]).")
	end
end
