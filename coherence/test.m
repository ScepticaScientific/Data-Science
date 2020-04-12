% This utility provides a test for the canonical coherence analysis 
% function 'getCanonicalCoherence()'.
clear all;

%% Inits
N = 3;              % Number of variates in the vector (multivariate) time series

fs = 1000.0;        % Discretisation frequency
t = [0.0 : 1.0 / fs : 1.0 - 1.0 / fs]';

fcommon = 200.0;    % Base frequency

% We generate a multivariate time series
ddx = repmat(1.0 * cos(2.0 * pi * t * fcommon) + 0.5 * cos(2.0 * pi * t * fcommon / 2.0) + 1.0 * cos(2.0 * pi * t * fcommon / 4.0), 1, N) + randn(length(t), N);

%% Computing
[evt, ev, freq] = getCanonicalCoherence(ddx);

%% Output
figure;
for i = 1 : N
    subplot(N, 1, i);
    semilogx(freq * fs, ev(:, i));
    if (i == N)
        xlabel('$\omega$');
    end
    ylabel(sprintf('$c_%d(\\omega)$', i));
    if (i == 1)
        title('Canonical Coherence Analysis (non-parametric estimation)');
    end
    set(gca, 'TickLabelInterpreter', 'latex');
end

figure;
semilogx(freq * fs, evt);
xlabel('$\omega$');
ylabel('$C(\omega)$');
title('Canonical Coherence Analysis (non-parametric estimation)');
set(gca, 'TickLabelInterpreter', 'latex');
