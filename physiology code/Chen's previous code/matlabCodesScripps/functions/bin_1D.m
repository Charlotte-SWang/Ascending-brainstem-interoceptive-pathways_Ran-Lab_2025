

% bin_1D
% bin original data to other bin size. only for 1D data(vector). 
% to bin 3D data ,see bin_3D
% yjzhu. 2008/03/28

function [new_data] = bin_1D(data,bin_size)


[dim_a dim_b] = size(data);

if (dim_a ~= 1) & (dim_b ~= 1) 
    Error('Only receive vector input!')
end;

bins = floor(length(data)./bin_size);

new_data = zeros(bins,1);

for i = 1:bins
    t_start = (i-1)*bin_size +1;
    t_end = t_start + bin_size -1;
    new_data(i) = mean(data(t_start:t_end));
end;
