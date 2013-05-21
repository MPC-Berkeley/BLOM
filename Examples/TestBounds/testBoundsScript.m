%% testOverLapping Bounds Case
horizon = 4;

%% test 1
set_param('testOverlappingBounds/Bound','lb','-5')
set_param('testOverlappingBounds/Bound','ub','100')
set_param('testOverlappingBounds/Bound1','lb','-3')
set_param('testOverlappingBounds/Bound1','ub','101')
[ModelSpec,block,stepVars,allVars] = BLOM_ExtractModel('testOverlappingBounds',horizon,1,1,1);
if isequal([stepVars.initLowerBound stepVars.interLowerBound stepVars.finalLowerBound],[-3 -3 -3])
    fprintf('The lower bound test passed\n')
else
    fprintf('The lower bound test failed\n')
end


if isequal([stepVars.initUpperBound stepVars.interUpperBound stepVars.finalUpperBound],[100 100 100])
    fprintf('The upper bound test passed\n')
else
    fprintf('The upper bound test failed\n')
end