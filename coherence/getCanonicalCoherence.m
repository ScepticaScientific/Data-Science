function [evt, ev, freq] = getCanonicalCoherence(xx)
% [evt, ev, freq] = getCanonicalCoherence(xx)
%
% This code implements the canonical coherence analysis of multivariate 
% data. The computation is performed using a non-parametric spectrum 
% estimation. For details, please refer to [1].
%
% At the input, 'xx' is a multivariate time series of second increments of 
% the physical observables 'x(t) = (x_1(t), ..., x_N(t))'.
%
% At the output, 'evt' is the total coherence coefficient over the 
% frequency range, 'ev' is the array of partial coherence coefficients over 
% the frequency range, 'freq' is the frequency range over which the 
% coherence spectrum is computed.
%
% REFERENCES:
% [1] A.A. Lyubushin, Data Analysis of Systems of Geophysical and 
%     Ecological Monitoring, Nauka, Moscow, 2007.
% [2] D.M. Filatov and A.A. Lyubushin, Physica A: Stat. Mech. Appl., 527 
%     (2019) 121309. DOI: 10.1016/j.physa.2019.121309.
%
% The end user is granted perpetual permission to reproduce, adapt, and/or 
% distribute this code, provided that an appropriate link is given to the 
% original repository it was downloaded from.

%% Auxiliaries
N = size(xx, 2);    % Number of variates in the vector time series

% Determining the frequency range
[~, freq] = cpsd(xx(:, 1), xx(:, 1), [], [], [], 1);
freq_len = length(freq);

%% Computing
ev = zeros(freq_len, N);
for ic = 1 : N
    %% We split the original N-variate signal into an (N - 1)-variate and a 
    % single-variate ones
    x = xx(:, [1 : ic - 1 ic + 1 : N]);     % (N - 1)-variate (i.e. vector) time series
    y = xx(:, ic);                          % Single-variate (i.e. scalar) time series

    %% We compute the spectral matrices (see [1, p. 111]) ...
    % ... for the first signal, ...
    Sxx = zeros(N - 1, N - 1, freq_len);
    for i = 1 : N - 1
        for j = 1 : N - 1
            Sxx(i, j, :) = cpsd(x(:, i), x(:, j), [], [], [], 1);
        end
    end

    % ... for the mixture of the first and second signals, ...
    Sxy = zeros(N - 1, 1, freq_len);
    for i = 1 : N - 1
        Sxy(i, 1, :) = cpsd(x(:, i), y(:, 1), [], [], [], 1);
    end
    
    Syx = zeros(1, N - 1, freq_len);
    for j = 1 : N - 1
        Syx(1, j, :) = cpsd(y(:, 1), x(:, j), [], [], [], 1);
    end
        
    % ... and for the second signal
    Syy = zeros(1, 1, freq_len);
    Syy(1, 1, :) = cpsd(y(:, 1), y(:, 1), [], [], [], 1);

    %% We compute the matrix whose eigenvalues are the measures of 
    % coherence and explicitly determine the maximum value of the coherence 
    % at each frequency
    for i = 1 : freq_len
        % We implement formula (2.3.1) from [1]
        U = inv(sqrtm(Sxx(:, :, i))) * Sxy(:, :, i) * inv(Syy(:, :, i)) * Syx(:, :, i) * inv(sqrtm(Sxx(:, :, i)));
        
        ev(i, ic) = max(real(eig(U)));
    end
end

% Finally we compute the total coherence frequency-wise
evt = prod(ev, 2);
