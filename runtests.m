clear all
import matlab.unittest.TestSuite

approxTests = TestSuite.fromPackage("approx");

suite = [approxTests];

disp(suite.run().table())
