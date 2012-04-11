function code = BLOM_FunctionCode(fcn)
% This function outputs the exception values used in BLOM to represent
% special transcendental (non-polynomial) functions. It takes one input
% argument, a string representing the desired function such as 'exp',
% 'log', etc. The output is a scalar double value that is used everywhere
% else in BLOM to represent these functions and trigger special behavior
% for function and gradient generation and evaluation.

switch fcn
   case 'exp'
      code = 1000000;
   case 'log'
      code = 1000001;
   otherwise
      % output a vector listing all the exception code values
      code = [1000000, 1000001];
end
