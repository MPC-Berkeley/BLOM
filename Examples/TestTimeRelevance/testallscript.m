% test script
simin = [0 1];

% testManyBound test script
horizon = 4;
[a,blockMB,stepVarsMB,allVarsMB] = BLOM_ExtractModel('testManyBound',horizon,1,1,1);

% testManyDelay test script
horizon=4;
[a,blockMD,stepVarsMD,allVarsMD] = BLOM_ExtractModel('testManyDelay',horizon,1,1,1);

% testSubsystem test script
horizon=4;
[a,blockS,stepVarsS,allVarsS] = BLOM_ExtractModel('testSubsystem',horizon,1,1,1);
