clear all
import matlab.unittest.TestSuite
tic

approxTests = TestSuite.fromPackage("approx");
folderTests = TestSuite.fromFolder("test");

suite = [approxTests, folderTests];
rng(1);
disp(suite.run().table())
toc
