% This utility provides tests for the canonical coherence analysis 
% functions 'getCanonicalCoherence()' and 'getCanonicalCoherenceW()'.
clear all;

%% Inits
testID = 2;

if (testID == 1)        % No coherence
    N = 3;              % Number of variates in the vector (multivariate) time series
    fs = 1000.0;        % Discretisation frequency
    t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';
    fcommon = 200.0;    % Base frequency
    
    ddx = zeros(length(t), N);
    ddx(:, 1) = cos(2.0 * pi * t * fcommon) + randn(length(t), 1);
    ddx(:, 2) = 0.5 * cos(2.0 * pi * t * fcommon / 2.0) + randn(length(t), 1);
    ddx(:, 3) = cos(2.0 * pi * t * fcommon / 4.0) + randn(length(t), 1);
elseif (testID == 2)    % Coherence at three different frequencies
    N = 3;              
    fs = 2000.0;       
    t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';
    fcommon = 200.0;    
    
    ddx = repmat(1.0 * cos(2.0 * pi * t * fcommon) + 0.5 * cos(2.0 * pi * t * fcommon / 2.0) + 1.0 * cos(2.0 * pi * t * fcommon / 4.0), 1, N) + randn(length(t), N);
elseif (testID == 3)    % Coherence at two different frequencies in time-shifted time ranges
    N = 3;
    fs = 1000.0;
    t = [0.0 : 1.0 / fs : 2.0 - 1.0 / fs]';
    fcommon1 = 10.0;
    fcommon2 = 50.0;
    
    ddx = zeros(length(t), N);
    ddx(:, 1) = cos(2.0 * pi * fcommon1 * t) .* (t >= 0.5 & t < 1.1) + ...
                cos(2.0 * pi * fcommon2 * t) .* (t >= 0.2 & t < 1.4) + 0.25 * randn(size(t));
    ddx(:, 2) = sin(2.0 * pi * fcommon1 * t) .* (t >= 0.6 & t < 1.2) + ...
                sin(2.0 * pi * fcommon2 * t) .* (t >= 0.4 & t < 1.6) + 0.35 * randn(size(t));
    ddx(:, 3) = sin(2.0 * pi * fcommon1 * t) .* (t >= 1.5 & t < 1.8) + ...
                sin(2.0 * pi * fcommon2 * t) .* (t >= 1.3 & t < 1.7) + 0.15 * randn(size(t));
end

%% Computing
[evt, ev, freq] = getCanonicalCoherence(ddx);
[evt_w, ev_w, freq_w] = getCanonicalCoherenceW(ddx);

%% Output
% Fourier-based CCA
figure;
for i = 1 : N
    subplot(N, 1, i);
    semilogx(freq * fs, ev(:, i));
    if (i == N)
        xlabel('\omega');
    end
    ylabel(sprintf('c_%d(\\omega)', i));
    if (i == 1)
        title('Canonical Coherence Analysis (non-parametric Fourier spectrum estimation)');
    end
end

figure;
semilogx(freq * fs, evt);
xlabel('\omega');
ylabel('C(\omega)');
title('Canonical Coherence Analysis (non-parametric Fourier spectrum estimation)');

% Wavelet-based CCA
figure;
for i = 1 : N
    subplot(N, 1, i);
    
    pcolor(t, freq_w * fs, ev_w(:, :, i));
    set(gca, 'YScale', 'log');
    shading flat;
    colorbar;
    if (i == N)
        xlabel('t');
    end
    ylabel('\omega');
    if (i == 1)
        title('Canonical Coherence Analysis (wavelet spectrum estimation)');
    end
    axis xy;
    axis tight;    
end

figure;
pcolor(t, freq_w * fs, evt_w);
set(gca, 'YScale', 'log');
shading flat;
colorbar;
xlabel('t');
ylabel('\omega');
title('Canonical Coherence Analysis (wavelet spectrum estimation)');
axis xy;
axis tight;
