function SolverStruct = BLOM_ExportToSolver(ModelSpec,solver,options)
%
% SolverStruct = BLOM_ExportToSolver(ModelSpec,solver,options)
%
% Prepares problem for export to a solver. Handles all tasks that are not
% data dependant. 
%
% Input:
%   ModelSpec -   Model structure generatated by BLOM_ExtractModel.
%   solver    -   Solver name {'fmincon','IPOPT')
%   options   -   options created by BLOM_optset function.
%


switch (solver)
    case 'fmincon'

        [cost_name costGrad_name eqconstr_name eqconstrGrad_name neqconstr_name neqconstrGrad_name ...
            Aeq Beq A B ] = CreateProblemCBs(ModelSpec.name,ModelSpec.all_names, ModelSpec.AAs ,  ModelSpec.Cs , ModelSpec.ineq,ModelSpec.cost,true);

        rehash; % needed to reload the newly generated functions to Matlab's memory
        
        pr.objective = @(x)GenericCostFun(x,str2func(cost_name),str2func(costGrad_name)) ;
        pr.nonlcon = @(x)GenericConstrFun(x,str2func(neqconstr_name),str2func(eqconstr_name),...
            str2func(neqconstrGrad_name),str2func(eqconstrGrad_name)) ;
        pr.Aineq = A;
        pr.bineq = B;
        pr.Aeq = Aeq;
        pr.beq = Beq;
        pr.lb = -inf(size(ModelSpec.A,2),1);
        pr.ub = inf(size(ModelSpec.A,2),1);
        pr.solver = 'fmincon';
        
        pr.options = optimset('GradObj','on','GradConstr','on','MaxIter',1000);


        SolverStruct.solver = solver;
        SolverStruct.pr  = pr;
        SolverStruct.name = ModelSpec.name;
        
%         idx = find(ModelSpec.ex_vars~=0 | ModelSpec.all_state_vars~=0 );
%         for i=1:length(idx)
% %             if (strcmp(mode,'fmincon'))
% %                 pr.Aeq(end+1,idx(i))=1;
% %                 pr.beq(end+1) = pr.x0(idx(i));
% %             end
%         end
                
    case 'IPOPT'
        idx = find(ModelSpec.ex_vars~=0 | ModelSpec.all_state_vars~=0 );
        %{
        k=0;
        for i=1:length(idx)
            k = k+1;
            fixed.AAs{k} = sparse(length(ModelSpec.all_names),2);
            fixed.AAs{k}(idx(i),1)=1;
            fixed.Cs{k}(1) = -1;
            fixed.Cs{k}(2) = nan;
        end
        fixed_old = fixed;
        %}
        fixed.AAs = repmat({spalloc(length(ModelSpec.all_names),2,1)}, 1, length(idx));
        for i=1:length(idx)
            fixed.AAs{i}(idx(i),1)=1;
        end
        fixed.Cs = repmat({[-1 nan]}, 1, length(idx));
        %if ~isequalwithequalnans(fixed, fixed_old)
        %    warning('mismatch in fixed struct')
        %end
        fixed.idx = idx;
        
        AAs = repmat({[]}, 1, length(ModelSpec.AAs));
        Cs = repmat({[]}, 1, length(ModelSpec.Cs));
        AAs{1} = ModelSpec.A(ModelSpec.eq_start_A:ModelSpec.eq_end_A, :);
        Cs{1} = ModelSpec.C(ModelSpec.eq_start_C:ModelSpec.eq_end_C, ...
            ModelSpec.eq_start_A:ModelSpec.eq_end_A);
        ineq.AAs = repmat({[]}, 1, length(ModelSpec.ineq.AAs));
        ineq.Cs = repmat({[]}, 1, length(ModelSpec.ineq.Cs));
        ineq.AAs{1} = ModelSpec.A(ModelSpec.ineq_start_A:ModelSpec.ineq_end_A, :);
        ineq.Cs{1} = ModelSpec.C(ModelSpec.ineq_start_C:ModelSpec.ineq_end_C, ...
            ModelSpec.ineq_start_A:ModelSpec.ineq_end_A);
        %{
        if ~isequal(vertcat(AAs{:}), vertcat(ModelSpec.AAs{:}))
            warning('mismatch in AAs')
        elseif ~isequal(blkdiag(Cs{:}), blkdiag(ModelSpec.Cs{:}))
            warning('mismatch in Cs')
        elseif ~isequal(vertcat(ineq.AAs{:}), vertcat(ModelSpec.ineq.AAs{:}))
            warning('mismatch in ineq.AAs')
        elseif ~isequal(blkdiag(ineq.Cs{:}), blkdiag(ModelSpec.ineq.Cs{:}))
            warning('mismatch in ineq.Cs')
        end
        %}
        
        CreateIpoptCPP(ModelSpec.name, ModelSpec.all_names, AAs, Cs, ...
            ineq, fixed, ModelSpec.cost);
        SolverStruct.solver = solver;
        SolverStruct.name = ModelSpec.name;
        SolverStruct.fixed_idx = idx;
    case 'linprog'
        SolverStruct.solver = solver;
        SolverStruct.name = ModelSpec.name;

        
        if BLOM_CheckIfLinear(ModelSpec.A) == false
            warning('The model is not linear, will be linearized around x=0');
        end
        Jac       = BLOM_EvalJacobian(ModelSpec.A, ModelSpec.C , zeros(size(ModelSpec.A,2),1));
        Constant  = BLOM_EvalPolyBlock(ModelSpec.A,ModelSpec.C, zeros(size(ModelSpec.A,2),1));
        SolverStruct.pr.f = Jac(1,:);
        SolverStruct.pr.Aineq  = Jac(ModelSpec.ineq_start_C:ModelSpec.ineq_end_C,:);
        SolverStruct.pr.bineq = -Constant(ModelSpec.ineq_start_C:ModelSpec.ineq_end_C);
        SolverStruct.pr.Aeq  = Jac(ModelSpec.eq_start_C:ModelSpec.eq_end_C,:);
        SolverStruct.pr.beq = -Constant(ModelSpec.eq_start_C:ModelSpec.eq_end_C);
        SolverStruct.pr.lb = -inf(size(ModelSpec.A,2),1);
        SolverStruct.pr.ub = inf(size(ModelSpec.A,2),1);
        SolverStruct.pr.solver = 'linprog';
        SolverStruct.pr.options = optimset;

        
        
end