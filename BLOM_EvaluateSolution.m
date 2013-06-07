function [Eq Ineq Cost] = BLOM_EvaluateSolution(ModelSpec,data)
%
% [Eq Ineq Cost] = BLOM_EvaluateSolution(ModelSpec,data)
%
% Evaluates the equality and inequality constraints and the cost value in
% the given point.
%
% Input:
%   ModelSpec -  Model data.
%   data      -  Solution vecotr either as a vector or as a structure.
%
% Output:
%   Eq        - Vector containing equality constraints values
%   Ineq      - Vector containing inequality constraints values
%   Cost      - Value of the cost function


if  isstruct(data)
    xvec = BLOM_ConvertStructToVector(ModelSpec.all_names,data);
else
    xvec = data;
end

Eq = [];
for i=1:length(ModelSpec.AAs)
    val = BLOM_EvalPolyBlock(ModelSpec.AAs{i}',ModelSpec.Cs{i},xvec);
    Eq = [Eq ; val ];
end

Ineq = [];
for i=1:length(ModelSpec.ineq.AAs)
    val = BLOM_EvalPolyBlock(ModelSpec.ineq.AAs{i}',ModelSpec.ineq.Cs{i},xvec);
    Ineq = [Ineq ; val ];
end

Cost =  BLOM_EvalPolyBlock(ModelSpec.cost.A ,ModelSpec.cost.C,xvec);


