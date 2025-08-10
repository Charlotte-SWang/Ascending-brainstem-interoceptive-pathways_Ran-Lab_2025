function [outputArg1] = calculateDistribution(inputdata,distributionBinSize,range)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nBin = (range(2) - range(1)) / distributionBinSize;
initialValue = range(1);
for iBin = 1:nBin
    distrib(iBin) = sum(inputdata > initialValue  + distributionBinSize * (iBin - 1) & inputdata <= initialValue + distributionBinSize * iBin);
end;
    
outputArg1 = distrib;
end

