function [evt, ev, freq] = getCanonicalCoherenceW(ddx, fs, isInfo)
% [evt, ev, freq] = getCanonicalCoherenceW(ddx, fs, isInfo)
%
% This code implements the canonical coherence analysis of multivariate 
% data. The computation is performed using a wavelet spectrum estimation. 
% For details, please refer to [1].
%
% At the input, 'ddx' is a multivariate time series of the second 
% increments of the physical observable 'ddx(t) = (ddx_1(t), ..., 
% ddx_N(t))', 'fs' is the sample rate (optional), 'isInfo' is the flag
% prescribing to output a message for each variate processed (optional).
%
% At the output, 'evt' is the total coherence spectrum (the coherence
% coefficient between 0.0 and 1.0), 'ev' is the array of the partial 
% coherence coefficients, while 'freq' is the frequency range over which 
% the coherence has been computed.
%
% REFERENCES:
% [1] A.A. Lyubushin, Data Analysis of Systems of Geophysical and 
%     Ecological Monitoring, Nauka, Moscow, 2007.
%
% The end user is granted perpetual permission to reproduce, adapt, and/or 
% distribute this code, provided that an appropriate link is given to the 
% original repository it was downloaded from.

    %% Auxiliaries
    if (nargin == 1)
        fs = 1.0;
        isInfo = false;
    elseif (nargin == 2)
        isInfo = false;
    end
    
    N = size(ddx, 2);    % Number of variates in the vector time series

    % Determining the frequency range ...
    [~, psd_freq] = pwelch(ddx(:, 1), [], [], [], fs);
    % ... for computing a consistent one for the CCWA
    [~, aux, freq] = wcoherence(ddx(:, 1), ddx(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
    freq_len = length(freq);
    
    t_len = size(aux, 2);

    %% Computing
    ev = zeros(freq_len, t_len, N);
    for ic = 1 : N
        %% We split the original N-variate signal into an (N - 1)-variate 
        % and a single-variate ones
        x = ddx(:, [1 : ic - 1 ic + 1 : N]);     % (N - 1)-variate (i.e. vector) time series
        y = ddx(:, ic);                          % Single-variate (i.e. scalar) time series

        %% We compute the spectral matrices ...
        % ... for the first signal, ...
        Sxx = zeros(N - 1, N - 1, freq_len, t_len);
        for i = 1 : N - 1
            for j = 1 : N - 1
                [~, Sxx(i, j, :, :), ~] = wcoherence(x(:, i), x(:, j), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
            end
        end

        % ... for the mixture of the first and second signals, ...
        Sxy = zeros(N - 1, 1, freq_len, t_len);
        for i = 1 : N - 1
            [~, Sxy(i, 1, :, :), ~] = wcoherence(x(:, i), y(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
        end

        Syx = zeros(1, N - 1, freq_len, t_len);
        for j = 1 : N - 1
            [~, Syx(1, j, :, :), ~] = wcoherence(y(:, 1), x(:, j), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
        end

        % ... and for the second signal
        Syy = zeros(1, 1, freq_len, t_len);
        [~, Syy(1, 1, :, :), ~] = wcoherence(y(:, 1), y(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);

        %% We compute the matrix, whose eigenvalues are PCA-related 
        % measures of coherence, and keep the maximum value
        %{
        % Way 1 (element-wise)
        for i = 1 : freq_len
            for j = 1 : t_len
                U = inv(sqrtm(Sxx(:, :, i, j))) * Sxy(:, :, i, j) * inv(Syy(:, :, i, j)) * Syx(:, :, i, j) * inv(sqrtm(Sxx(:, :, i, j)));

                ev(i, j, ic) = max(real(eig(U)));
            end
        end
        %}

        %
        % Way 2 (vectorised)
        Sxx_aux = toMatrixArray(Sxx);
        Sxx_aux = cellfun(@inv, Sxx_aux, 'UniformOutput', false);

        Sxy_aux = toMatrixArray(Sxy);

        Syy_aux = toMatrixArray(Syy);
        Syy_aux = cellfun(@inv, Syy_aux, 'UniformOutput', false);

        Syx_aux = toMatrixArray(Syx);

        aux = cellfun(@mtimes, Sxx_aux, Sxy_aux, 'UniformOutput', false);
        aux = cellfun(@mtimes, aux, Syy_aux, 'UniformOutput', false);
        U = cellfun(@mtimes, aux, Syx_aux, 'UniformOutput', false);
        
        aux = cellfun(@eig, U, 'UniformOutput', false);
        aux = cellfun(@real, aux, 'UniformOutput', false);
        ev(:, :, ic) = reshape(cellfun(@max, aux), freq_len, t_len);
        %}
        
        % Informational message to the standard output
        if (isInfo)
            fprintf('CCWA: variate %d of %d processed\n', ic, N);
        end
    end

    % Finally we compute the total coherence frequency-wise
    evt = prod(ev, 3) .^ (1.0 / N);
end

% This auxiliary function reshapes a 4-D array into a sequence of 2-D
% matrices
function res = toMatrixArray(mtrx)
    aux = reshape(mtrx, size(mtrx, 1), size(mtrx, 4) * size(mtrx, 3) * size(mtrx, 2)); 
    res = mat2cell(aux, size(mtrx, 1), ones(1, size(mtrx, 4) * size(mtrx, 3)) * size(mtrx, 2));
end    
