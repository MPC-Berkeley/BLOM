%> @file CombineFandG.m
%> @brief using user friendly P and K matrices for f(x) and g(x) in
%> f(x)/g(x) in General Polyblock, generates P and K matrix for the block
%>
%> @param P_f user friendly P matrix for f(x)
%> @param K_f user friendly K matrix for g(x)
%> @param P_g user friendly P matrix for f(x)
%> @param K_g user friendly K matrix for g(x)
%>
%> @retval P P matrix for the General Polyblock
%> @retval K K matrix for the General Polyblock
%========================================================================
function [P K] = CombineFandG(P_f,K_f,P_g,K_g)

P = [P_f, sparse(size(P_f,1),size(K_f,1)); ...
    repmat(P_g,size(K_f,1),1), kron(speye(size(K_f,1)),ones(size(P_g,1),1))];
if size(K_g,1) == 1
    K = [K_f, kron(speye(size(K_f,1)),-K_g)];
elseif size(K_g,1) == size(K_f,1)
    K_g_cell = num2cell(K_g, 2); % cell array of rows of K_g
    K = [K_f, -blkdiag(K_g_cell{:})];
else
    error('Number of rows in K matrices of f and g must be the same, or K_g must be a single row');
end
