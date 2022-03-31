% This utility provides tests for the canonical coherence analysis 
% functions 'getCanonicalCoherence()' and 'getCanonicalCoherenceW()'. Set 
% the testID parameter to run a certain test.
clear all;

%% Inits
testID = 4;

if (testID == 1)        % No coherence
    N = 3;              % Number of variates in the vector (multivariate) time series
    fs = 1000.0;        % Discretisation frequency
    t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';
    fcommon = 200.0;    % Base frequency
    
    ddx = zeros(length(t), N);
    ddx(:, 1) = cos(2.0 * pi * t * fcommon) + randn(length(t), 1);
    ddx(:, 2) = 0.5 * cos(2.0 * pi * t * fcommon / 2.0) + randn(length(t), 1);
    ddx(:, 3) = cos(2.0 * pi * t * fcommon / 4.0) + randn(length(t), 1);
elseif (testID == 2)    % Coherence at two different frequencies
    N = 2;              
    fs = 1000.0;       
    t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';
    fcommon = 200.0;    
    
    ddx = repmat(1.0 * cos(2.0 * pi * t * fcommon) + 1.0 * cos(2.0 * pi * t * fcommon / 4.0), 1, N) + randn(length(t), N);
elseif (testID == 3)    % Coherence at two different frequencies in time-shifted time ranges
    N = 2;
    fs = 1000.0;
    t = [0.0 : 1.0 / fs : 2.0 - 1.0 / fs]';
    fcommon1 = 10.0;
    fcommon2 = 50.0;
    
    ddx = zeros(length(t), N);
    ddx(:, 1) = cos(2.0 * pi * fcommon1 * t) .* (t >= 0.5 & t < 1.1) + ...
                cos(2.0 * pi * fcommon2 * t) .* (t >= 0.2 & t < 1.4) + 0.25 * randn(size(t));
    ddx(:, 2) = sin(2.0 * pi * fcommon1 * t) .* (t >= 0.6 & t < 1.2) + ...
                sin(2.0 * pi * fcommon2 * t) .* (t >= 0.4 & t < 1.6) + 0.35 * randn(size(t));
elseif (testID == 4)    % Coherence at one to two different frequencies in time-shifted time ranges
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

waveletSigma = 6.0;     % Default value is 6.0
energyLimit = 0.95;

%% Computing
[evt, ev, freq] = getCanonicalCoherence(ddx, fs, 1);
[evt_w, ev_w, freq_w, coi, timeBorders] = getCanonicalCoherenceW(ddx, fs, [t(floor(end / 4)) t(floor(4 * end / 5))], energyLimit, [], true);
[evt_wSgm, ev_wSgm, freq_wSgm, coiSgm, timeBordersSgm] = getCanonicalCoherenceW(ddx, fs, [t(floor(end / 4)) t(floor(4 * end / 5))], energyLimit, waveletSigma, true);

%% Output
% Fourier-based CCA
%
figure;
semilogx(freq, evt);
xlabel('f, Hz');
ylabel('C(f)');
title('Canonical Coherence Analysis (non-parametric Fourier spectrum estimation)');

figure;
for i = 1 : N
    subplot(N, 1, i);
    semilogx(freq, ev(:, i));
    if (i == N)
        xlabel('f, Hz');
    end
    ylabel(sprintf('c_%d(f)', i));
    if (i == 1)
        title('Canonical Coherence Analysis (non-parametric Fourier spectrum estimation)');
    end
end
%}

% Wavelet-based CCA, automatic time-frequency trade-off
%
figure;
pcolor(t, freq_w, evt_w);
set(gca, 'YScale', 'log');
shading flat;
colorbar;
%xlabel('t');
xlabel('t - b');
ylabel('f, Hz');
title('Wavelet CCA (automatic time-frequency trade-off), Total Coherence');
axis xy;
axis tight;
hold on;
plot(t, coi, 'w-.');

if (exist('timeBorders', 'var'))
    nb = size(timeBorders, 3);
    for ib = 1 : nb
        % We cut off the time moments (if any) which are out of the period 
        % of observation
        timeBorders(timeBorders(:, 1, ib) < t(1), 1, ib) = NaN;
        timeBorders(timeBorders(:, 2, ib) > t(end), 2, ib) = NaN;
        plot(timeBorders(:, 1, ib), freq_w, 'w--', timeBorders(:, 2, ib), freq_w, 'w--');
    end
end

figure;
for i = 1 : N
    subplot(N, 1, i);
    
    pcolor(t, freq_w, ev_w(:, :, i));
    set(gca, 'YScale', 'log');
    shading flat;
    colorbar;
    if (i == N)
        %xlabel('t');
        xlabel('t - b');
    end
    ylabel('f, Hz');
    if (i == 1)
        title('Wavelet CCA (automatic time-frequency trade-off), Partial Coherences');
    end
    axis xy;
    axis tight;
    hold on;
    plot(t, coi, 'w-.');
    
    if (exist('timeBorders', 'var'))
        for ib = 1 : nb
            plot(timeBorders(:, 1, ib), freq_w, 'w--', timeBorders(:, 2, ib), freq_w, 'w--');
        end
    end
end
%}

% Wavelet-based CCA, user-defined time-frequency trade-off
%
figure;
pcolor(t, freq_wSgm, evt_wSgm);
set(gca, 'YScale', 'log');
shading flat;
colorbar;
%xlabel('t');
xlabel('t - b');
ylabel('f, Hz');
title('Wavelet CCA (user-defined time-frequency trade-off), Total Coherence');
axis xy;
axis tight;
hold on;
plot(t, coiSgm, 'w-.');

if (exist('timeBordersSgm', 'var'))
    nb = size(timeBordersSgm, 3);
    for ib = 1 : nb
        % We cut off the time moments (if any) which are out of the period 
        % of observation
        timeBordersSgm(timeBordersSgm(:, 1, ib) < t(1), 1, ib) = NaN;
        timeBordersSgm(timeBordersSgm(:, 2, ib) > t(end), 2, ib) = NaN;
        plot(timeBordersSgm(:, 1, ib), freq_wSgm, 'w--', timeBordersSgm(:, 2, ib), freq_wSgm, 'w--');
    end
end

figure;
for i = 1 : N
    subplot(N, 1, i);
    
    pcolor(t, freq_wSgm, ev_wSgm(:, :, i));
    set(gca, 'YScale', 'log');
    shading flat;
    colorbar;
    if (i == N)
        %xlabel('t');
        xlabel('t - b');
    end
    ylabel('f, Hz');
    if (i == 1)
        title('Wavelet CCA (user-defined time-frequency trade-off), Partial Coherences');
    end
    axis xy;
    axis tight;
    hold on;
    plot(t, coiSgm, 'w-.');
    
    if (exist('timeBordersSgm', 'var'))
        for ib = 1 : nb
            plot(timeBordersSgm(:, 1, ib), freq_wSgm, 'w--', timeBordersSgm(:, 2, ib), freq_wSgm, 'w--');
        end
    end
end
%}