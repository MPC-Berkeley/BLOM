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

if (size(K_g,1)==1)
        P = [P_f zeros(size(P_f,1),size(K_f,1)) ; repmat(P_g,size(K_f,1),1) kron(eye( size(K_f,1)),ones(size(P_g,1),1))  ];
        K = [K_f kron(eye(size(K_f,1)),-K_g)];
elseif size(K_g,1)==size(K_f,1)
        P = [P_f zeros(size(P_f,1),size(K_f,1)) ; repmat(P_g,size(K_f,1),1) kron(eye( size(K_f,1)),ones(size(P_g,1),1))  ];
        K = [K_f zeros(size(K_g,1),size(K_g,2)*size(K_g,1))];
        for i=1:size(K_g,1)
            K(i,size(K_f,2)+1 + ((i-1)*size(K_g,2):(i*size(K_g,2)-1))) = -K_g(i,:);       
        end   
else
    error('number of lines in C matrices of f and g must be either the same, or K_f of g be a row matrix');
end
% input_type, output_type