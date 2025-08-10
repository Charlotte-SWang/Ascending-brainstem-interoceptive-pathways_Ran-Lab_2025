y = fft(a);
fs = sampleRate;
f = (0:length(y)-1)*fs/length(y);

plot(f,abs(y))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')
xlim([0,20])
ylim([0,10000])

plot(bandstop(a,[0.5 8],sampleRate))
hold on 
plot(a)