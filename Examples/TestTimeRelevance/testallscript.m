% test each model
testfailed = false;
testfailed = testfun('testManyBound')||testfailed;
testfailed = testfun('testManyDelay')||testfailed;
testfailed = testfun('testSubsystem')||testfailed;
testfailed = testfun('testSubSystemMIMO')||testfailed;
testfailed = testfun('testFromGoto')||testfailed;

if testfailed
    error('All tests run.  One or more tests failed')
end