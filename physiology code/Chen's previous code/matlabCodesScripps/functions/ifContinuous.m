function [outputArg1] = ifContinuous(A,threshold,nPoint)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
idx=A>=threshold;
ii1=strfind([0 idx 0],[0 1]);
ii2=strfind([0 idx 0],[1 0])-1;
ii=(ii2-ii1+1)>=nPoint;
outputArg1 = sum(ii);
end

