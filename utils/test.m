% This is a test utility for the computation of continuous wavelet power 
% spectra provided by the function 'getPowerSpectrumW()'.
clear all;

fs = 1000.0;                           
t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';
fcommon1 = 100.0;
fcommon2 = 10.0;
c1 = cos(2.0 * pi * t * fcommon1);
c2 = cos(2.0 * pi * t * fcommon2);

N = 2;
x = zeros(length(t), N);
x(:, 1) = c1 + randn(length(t), 1);
x(:, 2) = c2 + randn(length(t), 1);

%% Computing
[wPsx1x2, freqSx1x2, coix1x2] = getPowerSpectrumW(x, fs, 6.0);
[wPsx1x1, freqSx1x1, coix1x1] = getPowerSpectrumW(x(:, 1), fs, 6.0);

%% Output
figure;
pcolor(t, freqSx1x2, abs(wPsx1x2));
xlabel('Time, sec');
ylabel('Frequency, Hz');
shading interp;
set(gca, 'YScale', 'log');
hold on;
plot(t, coix1x2, 'w--');
title('Wavelet Power Cross-Spectrum');

figure;
pcolor(t, freqSx1x1, abs(wPsx1x1));
xlabel('Time, sec');
ylabel('Frequency, Hz');
shading interp;
set(gca, 'YScale', 'log');
hold on;
plot(t, coix1x1, 'w--');
title('Wavelet Power Spectrum');