clear all
init
import matlab.unittest.TestSuite
tic

approxTests = TestSuite.fromPackage("approx");
folderTests = TestSuite.fromFolder("test");

suite = [approxTests, folderTests];
rng(1);
results = suite.run();
disp(results.table())
toc

fprintf('\n')
results.assertSuccess();
disp('All tests passed!')
