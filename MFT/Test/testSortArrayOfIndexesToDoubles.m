% Class-based unit testing for SortArrayOfIndexesToDoubles
%
%   run these tests with two command-line commands:
%   - testCase = testSortArrayOfIndexesToDoubles;
%   - res = run(testCase);
%
classdef testSortArrayOfIndexesToDoubles < matlab.unittest.TestCase



methods (Test)
    function testSort (testCase)
       valueArr = [4.6, 38.5, 1e-3, 4.7, 6.8, 99.9, 1e5];
       indexArr = [6,3,4,2,7,1,5];
       minIx = 1;
       maxIx = 7;
       [actIndexArr]=SortArrayOfIndexesToDoubles(indexArr, valueArr, minIx, maxIx);
       [expIndexArr] = [3,1,4,5,2,6,7];
       testCase.verifyEqual(actIndexArr,expIndexArr)
    end
end

end