function CostFunction(block)
%MSFUNTMPL A template for an M-file S-function
%   The M-file S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl' with the name
%   of your S-function.  
%
%   It should be noted that the M-file S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more 
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%  
%   Copyright 2003-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  
  
%%
%% The setup method is used to setup the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.  
%%   
setup(block);
  
%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the S-function block's basic characteristics such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%% 
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

%    IneqSign = block.DialogPrm(1).Data;
%   C = block.DialogPrm(2).Data;
  
%    obj = gcbh;
%    set_param(obj,'UserData',IneqSign);
%   
  nIn =  1;
  nOut =  0;

%   Register number of ports
  block.NumInputPorts  = nIn;
  block.NumOutputPorts = nOut;
  
  % Setup port properties to be inherited or dynamic
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  % Override input port properties
  for i=1:nIn
      block.InputPort(i).DatatypeID  = 0;  % double
      block.InputPort(i).Complexity  = 'Real';
  end
  
  % Override output port properties
  for i=1:nOut
      block.OutputPort(1).DatatypeID  = 0; % double
      block.OutputPort(1).Complexity  = 'Real';
  end

  % Register parameters
  block.NumDialogPrms     = 0;
  block.DialogPrmsTunable = {};

  % Register sample times
  %  [0 offset]            : Continuous sample time
  %  [positive_num offset] : Discrete sample time
  %
  %  [-1, 0]               : Inherited sample time
  %  [-2, 0]               : Variable sample time
  block.SampleTimes = [0 0];
  
  %% -----------------------------------------------------------------
  %% Options
  %% -----------------------------------------------------------------
  % Specify if Accelerator should use TLC or call back into 
  % M-file
  block.SetAccelRunOnTLC(false);
  
  %% -----------------------------------------------------------------
  %% The M-file S-function uses an internal registry for all
  %% block methods. You should register all relevant methods
  %% (optional and required) as illustrated below. You may choose
  %% any suitable name for the methods and implement these methods
  %% as local functions within the same file.
  %% -----------------------------------------------------------------
    
  %% -----------------------------------------------------------------
  %% Register methods called during update diagram/compilation
  %% -----------------------------------------------------------------
  
  %% 
  %% CheckParameters:
  %%   Functionality    : Called in order to allow validation of
  %%                      block's dialog parameters. User is 
  %%                      responsible for calling this method
  %%                      explicitly at the start of the setup method
  %%   C-Mex counterpart: mdlCheckParameters
  %%
  block.RegBlockMethod('CheckParameters', @CheckPrms);

  %%
  %% SetInputPortSamplingMode:
  %%   Functionality    : Check and set input and output port 
  %%                      attributes specifying if port is operating 
  %%                      in sample-based or frame-based mode
  %%   C-Mex counterpart: mdlSetInputPortFrameData
  %%   (Signal Processing Blockset is required in order to set a port
  %%    to be frame-based)
  %%
  block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
  
  %%
  %% SetInputPortDimensions:
  %%   Functionality    : Check and set input and optionally output
  %%                      port dimensions
  %%   C-Mex counterpart: mdlSetInputPortDimensionInfo
  %%
  block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);

  %%
  %% SetOutputPortDimensions:
  %%   Functionality    : Check and set output and optionally input
  %%                      port dimensions
  %%   C-Mex counterpart: mdlSetOutputPortDimensionInfo
  %%
%   block.RegBlockMethod('SetOutputPortDimensions', @SetOutPortDims);
  
  %%
  %% SetInputPortDatatype:
  %%   Functionality    : Check and set input and optionally output
  %%                      port datatypes
  %%   C-Mex counterpart: mdlSetInputPortDataType
  %%
  block.RegBlockMethod('SetInputPortDataType', @SetInpPortDataType);
  
  %%
  %% SetOutputPortDatatype:
  %%   Functionality    : Check and set output and optionally input
  %%                      port datatypes
  %%   C-Mex counterpart: mdlSetOutputPortDataType
  %%
%   block.RegBlockMethod('SetOutputPortDataType', @SetOutPortDataType);
  
  %%
  %% SetInputPortComplexSignal:
  %%   Functionality    : Check and set input and optionally output
  %%                      port complexity attributes
  %%   C-Mex counterpart: mdlSetInputPortComplexSignal
  %%
  block.RegBlockMethod('SetInputPortComplexSignal', @SetInpPortComplexSig);
  
  %%
  %% SetOutputPortComplexSignal:
  %%   Functionality    : Check and set output and optionally input
  %%                      port complexity attributes
  %%   C-Mex counterpart: mdlSetOutputPortComplexSignal
  %%
%   block.RegBlockMethod('SetOutputPortComplexSignal', @SetOutPortComplexSig);
  
  %%
  %% PostPropagationSetup:
  %%   Functionality    : Setup work areas and state variables. Can
  %%                      also register run-time methods here
  %%   C-Mex counterpart: mdlSetWorkWidths
  %%
  block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);

  %% -----------------------------------------------------------------
  %% Register methods called at run-time
  %% -----------------------------------------------------------------
  
  %% 
  %% ProcessParameters:
  %%   Functionality    : Called in order to allow update of run-time
  %%                      parameters
  %%   C-Mex counterpart: mdlProcessParameters
  %%  
  block.RegBlockMethod('ProcessParameters', @ProcessPrms);

  %% 
  %% InitializeConditions:
  %%   Functionality    : Called in order to initialize state and work
  %%                      area values
  %%   C-Mex counterpart: mdlInitializeConditions
  %% 
  block.RegBlockMethod('InitializeConditions', @InitializeConditions);
  
  %% 
  %% Start:
  %%   Functionality    : Called in order to initialize state and work
  %%                      area values
  %%   C-Mex counterpart: mdlStart
  %%
  block.RegBlockMethod('Start', @Start);

  %% 
  %% Outputs:
  %%   Functionality    : Called to generate block outputs in
  %%                      simulation step
  %%   C-Mex counterpart: mdlOutputs
  %%
  block.RegBlockMethod('Outputs', @Outputs);

  %% 
  %% Update:
  %%   Functionality    : Called to update discrete states
  %%                      during simulation step
  %%   C-Mex counterpart: mdlUpdate
  %%
  block.RegBlockMethod('Update', @Update);

  %% 
  %% Derivatives:
  %%   Functionality    : Called to update derivatives of
  %%                      continuous states during simulation step
  %%   C-Mex counterpart: mdlDerivatives
  %%
  block.RegBlockMethod('Derivatives', @Derivatives);
  
  %% 
  %% Projection:
  %%   Functionality    : Called to update projections during 
  %%                      simulation step
  %%   C-Mex counterpart: mdlProjections
  %%
  block.RegBlockMethod('Projection', @Projection);
  
  %% 
  %% SimStatusChange:
  %%   Functionality    : Called when simulation goes to pause mode
  %%                      or continnues from pause mode
  %%   C-Mex counterpart: mdlSimStatusChange
  %%
  block.RegBlockMethod('SimStatusChange', @SimStatusChange);
  
  %% 
  %% Terminate:
  %%   Functionality    : Called at the end of simulation for cleanup
  %%   C-Mex counterpart: mdlTerminate
  %%
  block.RegBlockMethod('Terminate', @Terminate);

  %% -----------------------------------------------------------------
  %% Register methods called during code generation
  %% -----------------------------------------------------------------
  
  %%
  %% WriteRTW:
  %%   Functionality    : Write specific information to RTW file
  %%   C-Mex counterpart: mdlRTW
  %%
  block.RegBlockMethod('WriteRTW', @WriteRTW);
%endfunction

%% -------------------------------------------------------------------
%% The local functions below are provided for illustrative purposes
%% to show how you may implement the various block methods listed
%% above.
%% -------------------------------------------------------------------

function CheckPrms(block)
  
%   a = block.DialogPrm(1).Data;
%   if ~strcmp(class(a), 'double')
%     DAStudio.error('Simulink:block:invalidParameter');
%   end
%   
%endfunction

function ProcessPrms(block)

  block.AutoUpdateRuntimePrms;
 
%endfunction

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
%   block.OutputPort(1).SamplingMode  = fd;
  
%endfunction

function SetInpPortDims(block, idx, di)
  
  block.InputPort(idx).Dimensions = di;
%   block.OutputPort(1).Dimensions  = di;

%endfunction

function SetOutPortDims(block, idx, di)
  
%   block.OutputPort(idx).Dimensions = di;
  block.InputPort(1).Dimensions    = di;

%endfunction

function SetInpPortDataType(block, idx, dt)
  
  block.InputPort(idx).DataTypeID = dt;
%   block.OutputPort(1).DataTypeID  = dt;

%endfunction
  
function SetOutPortDataType(block, idx, dt)

%   block.OutputPort(idx).DataTypeID  = dt;
  block.InputPort(1).DataTypeID     = dt;

%endfunction  

function SetInpPortComplexSig(block, idx, c)
  
  block.InputPort(idx).Complexity = c;
%   block.OutputPort(1).Complexity  = c;

%endfunction 
  
function SetOutPortComplexSig(block, idx, c)

%   block.OutputPort(idx).Complexity = c;
  block.InputPort(1).Complexity    = c;

%endfunction 
    
function DoPostPropSetup(block)
  block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;

  %% Register all tunable parameters as runtime parameters.
  block.AutoRegRuntimePrms;

%endfunction

function InitializeConditions(block)
%endfunction

function Start(block)

  block.Dwork(1).Data = 0;
   
%endfunction

function WriteRTW(block)
  
   block.WriteRTWParam('matrix', 'M',    [1 2; 3 4]);
   block.WriteRTWParam('string', 'Mode', 'Auto');
   
%endfunction

function Outputs(block)

%   A = block.DialogPrm(1).Data;
%   C = block.DialogPrm(2).Data;
%   nIn =  size(A,2)-1;
%   nOut =  size(C,1);
% 
% 
%   
%   for i=1:nIn
%       X(i) = block.InputPort(1).Data; 
%   end
%   
%   y = 0;
%   for j=1:size(A,1)-1
%       Pr = C(j);
%       for i=1:nIn
%          Pr = Pr*X(i)^A(i,j);
%       end
%       y = y + Pr;
%   end
%   block.OutputPort(1).Data = y;
  
%endfunction

function Update(block)
  
  block.Dwork(1).Data = block.InputPort(1).Data;
  
%endfunction

function Derivatives(block)
%endfunction

function Projection(block)
%endfunction

function SimStatusChange(block, s)
  
  if s == 0
    disp('Pause has been called');
  elseif s == 1
    disp('Continue has been called');
  end
  
%endfunction
    
function Terminate(block)
%endfunction
 
