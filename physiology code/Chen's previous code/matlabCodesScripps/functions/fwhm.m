

function width = fwhm(x,y,lev50)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)

% y = furtherSmoothDeConvRes{1}(resp_zone{1},1);
% x = resp_zone{1};
y = y / max(y);
N = length(y);
% lev50 = 0.5;
% if y(1) < lev50                  % find index of center (max or min) of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
    disp('Pulse Polarity = Positive')
% else
%     [garbage,centerindex]=min(y);
%     Pol = -1;
%     disp('Pulse Polarity = Negative')
% end
% i = 2;
% while sign(y(i)-lev50) == sign(y(i-1)-lev50)
%     i = i+1;
% end                                   %first crossing is between v(i-1) & v(i)
for i = 2:N
    if sign(y(i)-lev50) ~= sign(y(i-1)-lev50)
        break;
    end;
end;
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
tleadx = i;
i = centerindex+1;                    %start search for next crossing at center
% while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))   
%     i = i+1;
% end
for j = i+1:N
    if sign(y(j)-lev50) ~= sign(y(j-1)-lev50)
        break;
    end;
end;
ttrailx = j;
i = j;  
if i ~= N
    Ptype = 1;  
    disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2; 
    disp('Step-Like Pulse, no second edge')
    ttrail = N;
    width = x(end) - tlead;
end
% figure;
% plot(y);
% hold on;
% plot(tleadx, y(tleadx), 'ro');
% plot(ttrailx, y(ttrailx), 'ro');
% 
% 
% 
