warning off;
%% initialize
simin = [0 1];
horizon = 1;
testFailed = false;
failedTests = '';

%% compile
testConvert2Polyblock([],[],[],'compile');

%% test sin
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function7');
Pdes = sparse([1 2 3 4], [1 2 3 4], [BLOM_FunctionCode('sin'), BLOM_FunctionCode('sin'), 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Sin'];
end

%% test cos
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function8');
Pdes = sparse([1 2 3 4], [1 2 3 4], [BLOM_FunctionCode('cos'), BLOM_FunctionCode('cos'), 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Cos'];
end

%% test tan
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function9');
Pdes = sparse([1 2 3 4 3 4], [1 2 3 4 1 2], [BLOM_FunctionCode('sin'), BLOM_FunctionCode('sin'), 1, 1, BLOM_FunctionCode('cos'), BLOM_FunctionCode('cos')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Tan'];
end

%% test sincos
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function10');
Pdes = sparse(1:8, 1:8, [BLOM_FunctionCode('sin'), BLOM_FunctionCode('sin'), BLOM_FunctionCode('cos'), BLOM_FunctionCode('cos') 1 1 1 1]);
Kdes = [speye(4) -speye(4)];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', SinCos'];
end

%% test asin
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function11');
Pdes = sparse([1 2 3 4], [1 2 3 4], [1 1 BLOM_FunctionCode('sin'), BLOM_FunctionCode('sin')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Asin'];
end

%% test acos
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function12');
Pdes = sparse([1 2 3 4], [1 2 3 4], [1 1 BLOM_FunctionCode('cos'), BLOM_FunctionCode('cos')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Acos'];
end

%% test atan
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Trigonometric Function13');
Pdes = sparse([1 2 3 4 1 2], [1 2 3 4 3 4], [1 1 BLOM_FunctionCode('sin'), BLOM_FunctionCode('sin'),BLOM_FunctionCode('cos'), BLOM_FunctionCode('cos')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Atan'];
end

%% test hypot
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function20');
Pdes = 2*speye(6);
Kdes = [1 0 1 0 -1 0; 0 1 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Hypot'];
end

%% test unary minus
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Unary Minus1');
Pdes = speye(4);
Kdes = [1 0 1 0; 0 1 0 1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Unary Minus'];
end

%% test square
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function11');
Pdes = sparse(1:4, 1:4, [2 2 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Square'];
end

%% test sqrt
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function14');
Pdes = sparse(1:4, 1:4, [.5 .5 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Sqrt'];
end

%% test transpose 
% [P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function21');
% Pdes = sparse(1:4, 1:4, [1 1 1 1]);
% Kdes = [1 0 -1 0; 0 1 0 -1];
% if ~isequal(P, Pdes) || ~isequal(K, Kdes)
%    testFailed = true;
%    failedTests = [failedTests ', Transpose'];
% end

%% test reciprocal
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function12');
Pdes = sparse(1:4, 1:4, [-1 -1 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Reciprocal'];
end

%% test natural log
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function15');
Pdes = sparse([1 2 3 4], [1 2 3 4], [1 1 BLOM_FunctionCode('exp'), BLOM_FunctionCode('exp')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Natural Log'];
end

%% test magnitude^2
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function16');
Pdes = sparse(1:4, 1:4, [2 2 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Magnitude^2'];
end

%% test gain
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Gain1');
Pdes = speye(4);
Kdes = [3 0 -1 0; 0 3 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Gain'];
end

%% test 10^u  
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function19');
Pdes = sparse(1:4, 1:4, [1 1 BLOM_FunctionCode('log'), BLOM_FunctionCode('log')]);
Kdes = [log(10) 0 -1 0; 0 log(10) 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', 10^u'];
end

%% test 1/sqrt 
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function17');
Pdes = sparse(1:4, 1:4, [-.5 -.5 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', 1/Sqrt'];
end

%% test log10  
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Math Function18');
Pdes = sparse(1:4, 1:4, [1/log(10) 1/log(10) BLOM_FunctionCode('exp'), BLOM_FunctionCode('exp')]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', log10'];
end

%% test bias 
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Bias1');
Pdes = sparse(1:4, 1:4, [1 1 1 1], 5, 4);
Kdes = [1 0 -1 0 1; 0 1 0 -1 1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Bias'];
end

%% test add/sub 
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Add1');
Pdes = speye(10);
Kdes = [1 0 1 0 -1 0 -1 0 -1 0; 0 1 0 1 0 -1 0 -1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Add/Sub'];
end


%% test sum 
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Sum1');
Pdes = speye(10);
Kdes = [1 0 1 0 -1 0 -1 0 -1 0; 0 1 0 1 0 -1 0 -1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Sum'];
end

%% test product/divide
[P K] = BLOM_Convert2Polyblock('testConvert2Polyblock/Product1');
Pdes = sparse([1 1 1 1 2 2 2 2 3 4], [1 3 5 7 2 4 6 8 9 10], [1 1 -1 -1 1 1 -1 -1 1 1]);
Kdes = [1 0 -1 0; 0 1 0 -1];
if ~isequal(P, Pdes) || ~isequal(K, Kdes)
   testFailed = true;
   failedTests = [failedTests ', Product'];
end
%% terminate

if testFailed
    disp(['Following tests failed: ' failedTests(3:end)])
else
    disp('All Tests Pass')
end

testConvert2Polyblock([],[],[],'term');

warning on;