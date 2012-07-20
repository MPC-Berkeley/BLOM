function [new old] = BLOM_UpdateOldIpoptDat(folder)
% This function converts old-formulation Ipopt data files in the input
% folder to new-formulation (all constraints as bounds) data files for
% the TAO interface, and probably future versions of the Ipopt interface.

% If called with no input arguments, use folder = pwd
if nargin == 0
    folder = pwd;
end

% Begin by loading all old-format data files
skip_Fixed_and_X0 = false;
for datfiles = {'params', 'testFixed', 'testX0'}
    if exist([folder filesep datfiles{1} '.dat'], 'file')
        % in case Fixed and X0 not generated for this dataset
        old.(datfiles{1}) = load([folder filesep datfiles{1} '.dat']);
    else
        skip_Fixed_and_X0 = true;
    end
end
old.n = old.params(1); % number of optimization variables in old format
old.n_fixed = old.params(2); % number of 'fixed' variables in old format
old.n_constr = old.params(3); % number of constraints in old format
old.nnz_jac_g = old.params(4); % nonzeros in Jacobian in old format
old.nnz_h_lag = old.params(5); % nonzeros in Lagrangian Hessian in old format
old.n_ineq = old.params(6); % number of inequalities in old format
for txtfiles = {'A', 'C', 'FixedStruct', 'HessianStruct', 'JacobianStruct', 'LambdaStruct'}
    old.(txtfiles{1}) = BLOM_LoadDataFile([folder filesep txtfiles{1} '.txt'], 'tripletmat_ascii');
end
old.A = old.A(any(old.C,1),:); % remove unused terms
old.C = old.C(:,any(old.C,1));
% Bounds on optimization variables
old.lb = -inf(old.n, 1); % preallocate lower bound vector
old.ub = inf(old.n, 1); % preallocate upper bound vector
new.FixedVars = (old.FixedStruct' ~= 0); % logical column vector
if ~skip_Fixed_and_X0
    old.lb(new.FixedVars) = old.testFixed; % lb = ub = fixed
    old.ub(new.FixedVars) = old.testFixed;
    old.ineq0 = BLOM_EvalPolyBlock(old.A, old.C(2:old.n_ineq+1,:), old.testX0);
    new.X0 = [old.testX0; -old.ineq0]; % old format inequality constraints
    % are nonpositive, new format slack variables are nonnegative
end

% New format includes only equality constraints, with lower and upper
% bounds on optimization variables. Have to introduce a slack variable
% for each inequality constraint from the old formulation.
new.n_vars = old.n + old.n_ineq; % number of optimization variables in new format
new.n_fixed = old.n_fixed; % number of 'fixed' variables in new format
new.n_constr = old.n_constr; % number of constraints in new format
new.FixedVars = [new.FixedVars; false(old.n_ineq, 1)]; % slack variables aren't fixed
new.lb = [old.lb; zeros(old.n_ineq, 1)]; % slack variables are nonnegative
new.ub = [old.ub; inf(old.n_ineq, 1)]; % slack variables are nonnegative
new.P = blkdiag(old.A, speye(old.n_ineq)); % exponent matrix
new.K = [old.C, [sparse(1, old.n_ineq); ...
    speye(old.n_constr, old.n_ineq)]]; % coefficient matrix
% Every old inequality constraint is now an equality constraint:
% old inequality value + new slack variable = 0
new.n_terms = size(new.P, 1); % number of polynomial terms in new format
new.n_lb = nnz(isfinite(new.lb)); % number of finite lower bounds in new format
new.n_ub = nnz(isfinite(new.ub)); % number of finite upper bounds in new format
if skip_Fixed_and_X0
    % n_lb and n_ub should include n_fixed, but if there was not an
    % existing Fixed.dat file to get the values from, the bound values
    % for the FixedVars will not have been set yet
    new.n_lb = new.n_lb + new.n_fixed;
    new.n_ub = new.n_ub + new.n_fixed;
end

% Now save the new format to data files
params_file = fopen([folder filesep 'params.txt'], 'w+');
fprintf(params_file, [ ...
    'Number.of.optimization.variables...n_vars: %d \n' ...
    'Number.of.constraints............n_constr: %d \n' ...
    'Number.of.polynomial.terms........n_terms: %d \n' ...
    'Number.of.finite.lower.bounds........n_lb: %d \n' ...
    'Number.of.finite.upper.bounds........n_ub: %d \n' ...
    'Number.of.fixed.variables.........n_fixed: %d \n' ], ...
    new.n_vars, new.n_constr, new.n_terms, new.n_lb, new.n_ub, new.n_fixed);
frewind(params_file)
new.params = fscanf(params_file, '%*s %d');
fclose(params_file);
% Save exponent matrix to P.tripletMat
BLOM_SaveDataFile(new.P, [folder filesep 'P.tripletMat'], 'tripletmat_binary');
%BLOM_SaveDataFile(new.P, [folder filesep 'P.txt'], 'tripletmat_ascii');
% Save coefficient matrix to K.tripletMat
BLOM_SaveDataFile(new.K, [folder filesep 'K.tripletMat'], 'tripletmat_binary');
%BLOM_SaveDataFile(new.K, [folder filesep 'K.txt'], 'tripletmat_ascii');
% Save bounds that are not 'fixed' (External/IC) variables in lb.sparseVec and ub.sparseVec
lb_nonfixed = new.lb;
lb_nonfixed(new.FixedVars) = -inf;
BLOM_SaveDataFile(lb_nonfixed, [folder filesep 'lb.sparseVec'], 'sparsevec_binary', -inf);
%BLOM_SaveDataFile(lb_nonfixed, [folder filesep 'lb.txt'], 'sparsevec_ascii', -inf);
ub_nonfixed = new.ub;
ub_nonfixed(new.FixedVars) = inf;
BLOM_SaveDataFile(ub_nonfixed, [folder filesep 'ub.sparseVec'], 'sparsevec_binary', inf);
%BLOM_SaveDataFile(ub_nonfixed, [folder filesep 'ub.txt'], 'sparsevec_ascii', inf);
if ~skip_Fixed_and_X0
    % Save 'fixed' (External/IC) variables in fixed.sparseVec
    fixed = nan(new.n_vars, 1);
    fixed(new.FixedVars) = old.testFixed;
    BLOM_SaveDataFile(fixed, [folder filesep 'fixed.sparseVec'], 'sparsevec_binary', nan);
    %BLOM_SaveDataFile(fixed, [folder filesep 'fixed.txt'], 'sparsevec_ascii', nan);
    % Save optimization variable guess in X0.denseVec
    BLOM_SaveDataFile(new.X0, [folder filesep 'X0.denseVec'], 'densevec_binary');
    %BLOM_SaveDataFile(new.X0, [folder filesep 'X0.txt'], 'densevec_ascii');
end
