function varargout = getCanonicalCoherenceW(ddx, fs, timesOfInterest, energyThreshold, waveletSigma, isInfo)
% varargout = getCanonicalCoherenceW(ddx, fs, timesOfInterest, energyThreshold, waveletSigma, isInfo)
%
% This code implements canonical coherence analysis of multivariate 
% data. The computation is performed using a wavelet spectrum estimation. 
% For details, please refer to [1-4].
%
% At the input:
%   - 'ddx' is a multivariate stationary time series (usually the second 
%     derivative of a physical observable 'x(t)') 'ddx(t) = (ddx_1(t), ..., 
%     ddx_N(t))'; the variates are provided column-wise
%   - 'fs' is the sampling rate (optional)
%   - 'waveletSigma' is the time-frequency resolution (optional)
%   - 'isInfo' is the flag prescribing to output a message for each variate 
%     processed (optional).
% 
% Besides, for an array of time moments this function allows to determine 
% the left and right borders, out of which the daughter wavelet, while 
% moving in time at a fixed scale, does not involve those time moments' 
% data samples for the computation of CWT:
%   - 'timesOfInterest' is a row vector of the time moments (these are 
%     usually extreme events reflected in the data) for which the left and 
%     right borders are wanted (optional)
%   - 'energyThreshold' is the probability with which the daughter wavelet
%     does not affect the time moments' data samples (optional).
%
% At the output:
%   - 'evt' is the total coherence spectrum (the coherence coefficient 
%     between 0.0 and 1.0)
%   - 'ev' is the array of the partial coherence coefficients
%   - 'freq' is the scale-related frequency range over which the coherence 
%     has been computed
%   - 'coi' is the cone of influence
%   - 'timeBorders' are borders of the energy cones out of which the above-
%     specified time moments' data samples are not affected by the daughter 
%     wavelet (optional).
%
% REFERENCES:
% [1] A.A. Lyubushin, Data Analysis of Systems of Geophysical and 
%     Ecological Monitoring, Nauka, Moscow, 2007.
% [2] J. Ashmead, Quanta, 1 (2012) 58-70.
% [3] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998)
%     61-78.
% [4] Michael X. Cohen, Parameters of Morlet wavelet (time-frequency
%     trade-off), https://www.youtube.com/watch?v=LMqTM7EYlqY
%
% The end user is granted perpetual permission to reproduce, adapt, and/or 
% distribute this code, provided that an appropriate link is given to the 
% original repository it was downloaded from.

    %% Auxiliaries
    if (nargin == 1)
        fs = 1.0;
        timesOfInterest = [];
        energyThreshold = [];
        waveletSigma = [];
        isInfo = false;        
    elseif (nargin == 2)
        timesOfInterest = [];
        energyThreshold = [];
        waveletSigma = [];
        isInfo = false;        
    elseif (nargin == 3)
        energyThreshold = [];
        waveletSigma = [];
        isInfo = false;        
    elseif (nargin == 4)
        waveletSigma = [];
        isInfo = false;
    elseif (nargin == 5)
        isInfo = false;
    end

    N = size(ddx, 2);    % Number of variates in the vector time series

    % Determining the frequency range for subsequent computing a 
    % consistent one for the CCWA
    [~, psd_freq] = pwelch(ddx(:, 1), [], [], [], fs);            
    
    % We obtain the cone of influence ...
    if (isempty(waveletSigma))
        [~, aux, freq, coi] = wcoherence(ddx(:, 1), ddx(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);     % Automatic time-frequency trade-off
    else
        [aux, freq, coi] = getPowerSpectrumW(ddx(:, [1 1]), fs, waveletSigma, [psd_freq(1) psd_freq(end)]);             % User-defined time-frequency trade-off
    end
    freq_len = length(freq);
    % ... and correct the values outside the minimum frequency
    coi(coi < freq(end)) = freq(end);

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
        if (isempty(waveletSigma))
            % Automatic time-frequency trade-off
            for i = 1 : N - 1
                for j = 1 : N - 1
                    [~, Sxx(i, j, :, :), ~] = wcoherence(x(:, i), x(:, j), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
                end
            end
        else
            % User-defined time-frequency trade-off
            for i = 1 : N - 1
                for j = 1 : N - 1
                    [Sxx(i, j, :, :), ~, ~] = getPowerSpectrumW(x(:, [i j]), fs, waveletSigma, [psd_freq(1) psd_freq(end)]);
                end
            end
        end

        % ... for the mixture of the first and second signals, ...
        Sxy = zeros(N - 1, 1, freq_len, t_len);
        if (isempty(waveletSigma))
            % Automatic time-frequency trade-off
            for i = 1 : N - 1
                [~, Sxy(i, 1, :, :), ~] = wcoherence(x(:, i), y(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
            end
        else
            % User-defined time-frequency trade-off
            for i = 1 : N - 1
                [Sxy(i, 1, :, :), ~, ~] = getPowerSpectrumW([x(:, i) y(:, 1)], fs, waveletSigma, [psd_freq(1) psd_freq(end)]);
            end
        end

        Syx = zeros(1, N - 1, freq_len, t_len);
        if (isempty(waveletSigma))
            % Automatic time-frequency trade-off
            for j = 1 : N - 1
                [~, Syx(1, j, :, :), ~] = wcoherence(y(:, 1), x(:, j), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
            end
        else
            % User-defined time-frequency trade-off
            for j = 1 : N - 1
                [Syx(1, j, :, :), ~, ~] = getPowerSpectrumW([y(:, 1) x(:, j)], fs, waveletSigma, [psd_freq(1) psd_freq(end)]);
            end
        end

        % ... and for the second signal
        Syy = zeros(1, 1, freq_len, t_len);
        if (isempty(waveletSigma))
            % Automatic time-frequency trade-off
            [~, Syy(1, 1, :, :), ~] = wcoherence(y(:, 1), y(:, 1), fs, 'FrequencyLimits', [psd_freq(1) psd_freq(end)]);
        else
            % User-defined time-frequency trade-off
            [Syy(1, 1, :, :), ~, ~] = getPowerSpectrumW([y(:, 1) y(:, 1)], fs, waveletSigma, [psd_freq(1) psd_freq(end)]);
        end

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
            if (isempty(waveletSigma))
                fprintf('CCWA (automatic time-frequency trade-off): variate %d of %d processed\n', ic, N);
            else
                fprintf('CCWA (user-defined time-frequency trade-off): variate %d of %d processed\n', ic, N);
            end
        end
    end

    % Finally we compute the total coherence frequency-wise
    evt = prod(ev, 3) .^ (1.0 / N);
    
    varargout{1} = evt;
    varargout{2} = ev;
    varargout{3} = freq;
    varargout{4} = coi;

    % For the particular time moments (if any) we obtain the energy cones
    if (~isempty(energyThreshold))
        timeBorders = getBorders(freq, energyThreshold, timesOfInterest, fs, waveletSigma);
        varargout{5} = timeBorders;
    end
end

% This auxiliary function reshapes a 4-D array into a sequence of 2-D
% matrices
function res = toMatrixArray(mtrx)
    aux = reshape(mtrx, size(mtrx, 1), size(mtrx, 4) * size(mtrx, 3) * size(mtrx, 2)); 
    res = mat2cell(aux, size(mtrx, 1), ones(1, size(mtrx, 4) * size(mtrx, 3)) * size(mtrx, 2));
end    

% This function computes the time borders (or energy cones) out of which 
% the moving daughter wavelet does not affect the specific time moments
function timeBorders = getBorders(freq, energyThreshold, timesOfInterest, fs, waveletSigma)
    if (isempty(waveletSigma))
        waveletSigma = 6.0;                     % The Morlet wavelet parameter 'sigma'
    end

    waveletFreq = waveletSigma / (2.0 * pi);    % The frequency of the sine harmonic of the Morlet wavelet, a constant, in hertz (or in cycles per second)    

    scales = waveletFreq * fs ./ freq * 1.0;    % In relative units
    t = 1.0 * erfinv(energyThreshold);          % In samples (samples * relative units)
    t_deltas = t * scales / fs;                 % In seconds

    N = length(timesOfInterest);
    timeBorders = zeros(length(freq), 2, N);
    for i = 1 : N
        timeBorders(:, :, i) = [timesOfInterest(i) - t_deltas timesOfInterest(i) + t_deltas];
    end
end
