function code = BLOM_FunctionCode(fcn)
% This function outputs the exception values used in BLOM to represent
% special transcendental (non-polynomial) functions. It takes one input
% argument, a string representing the desired function such as 'exp',
% 'log', etc. The output is a scalar double value that is used everywhere
% else in BLOM to represent these functions and trigger special behavior
% for function and gradient generation and evaluation. If no inputs are
% given, outputs a vector of all recognized exception code values.

codes_struct.exp  = 1e20;
codes_struct.log  = 2e20;
codes_struct.sin  = 3e20;
codes_struct.cos  = 4e20;
codes_struct.tanh = 5e20;
codes_struct.atan = 6e20;

if nargin == 0
    % return vector of all code values
    code = cellfun(@(s) codes_struct.(s), fieldnames(codes_struct));
else
    if isfield(codes_struct, fcn)
        code = codes_struct.(fcn);
    else
        if isempty(fcn) || strcmpi(fcn,'all_codes')
            % return vector of all code values
            code = cellfun(@(s) codes_struct.(s), fieldnames(codes_struct));
        elseif strcmpi(fcn,'codes_struct')
            % return structure, including new field for vector of all codes
            codes_struct.all_codes = cellfun(@(s) codes_struct.(s), ...
                fieldnames(codes_struct));
            code = codes_struct;
        else
            error(['Function ' fcn ' not recognized'])
        end
    end
end


% notes to self:
% hex2num('43d74910d52d3051') is the largest double for which (1-eps/2)^x > 0
% hex2num('43c62e42fefa39ef') is the largest double for which (1+eps)^x < inf
% hex2num('c3d62e42fefa39ee') is the smallest double for which (1-eps/2)^x < inf
% hex2num('c3c74910d52d3052') is the smallest double for which (1+eps)^x > 0

