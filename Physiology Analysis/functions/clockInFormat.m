function [outputArg1] = clockInFormat
%outputt clock in sring format yyyy-mm-dd hh-mm
junk1 = clock;
for i = 1:5
    if junk1(i) < 10
        junk{i} = ['0' num2str(junk1(i))];
    else
        junk{i} = num2str(junk1(i));
    end
end
outputArg1 = sprintf('%s-%s-%s %s-%s',junk{1},junk{2},junk{3},junk{4},junk{5});
end

