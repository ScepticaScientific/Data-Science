function [Ps, freqS, coi] = getPowerSpectrumW(x, fs, waveletSigma, preFreqs)
% [Ps, freqS, coi] = getPowerSpectrumW(x, fs, waveletSigma, preFreqs)
%
% This code implements the computation of the continuous wavelet auto- or 
% cross- power spectrum. It is based on MATLAB's function 'cwtft()' (see 
% also 'wcoherence()' there). A Morlet wavelet is employed. For details, 
% please refer to [1-3].
%
% At the input:
%   - 'x' is a uni- or bivariate time series of the physical observables 
%     'x_1(t)' and, optionally, 'x_2(t)'
%   - 'fs' is the sampling rate (optional)
%   - 'waveletSigma' is the Morlet wavelet's time-frequency resolution 
%     (optional)
%   - 'preFreqs' are the preset scale-related frequency limits to perform 
%     CWT within (optional; it will be computed automatically if not 
%     provided).
%
% If 'x_2(t)' is omitted then the auto-spectrum of 'x_1(t)' is computed; 
% otherwise the cross-spectrum between 'x_1(t)' and 'x_2(t)' is computed. 
% The time series are to be provided column-wise for each variate.
%
% At the output:
%   - 'Ps' is the power spectrum
%   - 'freqS' is the scale-related frequency range over which the spectrum 
%     is computed
%   - 'coi' is the cone of influence.
%
% REFERENCES:
% [1] J. Ashmead, Quanta, 1 (2012) 58-70.
% [2] C. Torrence and G.P. Compo, Bull. Am. Meteorol. Soc., 79 (1998)
%     61-78.
% [3] M.X. Cohen, Parameters of Morlet wavelet (time-frequency trade-off), 
%     https://www.youtube.com/watch?v=LMqTM7EYlqY
%
% The end user is granted perpetual permission to reproduce, adapt, and/or
% distribute this code, provided that an appropriate link is given to the
% original repository it was downloaded from.

    %% Auxiliaries
    if (nargin == 1)
        fs = 1.0;
        waveletSigma = 6.0;
        preFreqs = [];
    elseif (nargin == 2)
        waveletSigma = 6.0;
        preFreqs = [];
    elseif (nargin == 3)
        preFreqs = [];
    end

    % We precompute the scale-related frequency limits to subsequently 
    % compute the CWT within ...
    if (isempty(preFreqs))
        [~, preFreqs] = pwelch(x(:, 1), [], [], [], fs);
    end
    % ... and prepare a padding to later add at both sides of the signals
    x_len_adjusted = adjustSignalLength(size(x, 1));

    % We compute the continuous Fourier transform of the Morlet wavelet
    [wft, freqS, scales, coi] = MorletCWT(size(x, 1), fs, [preFreqs(1) preFreqs(end)], x_len_adjusted, waveletSigma);

    % Now we compute the CWT of the first time series ...
    cwtX = getCWT(x(:, 1), x_len_adjusted, wft);
    
    if (size(x, 2) == 1)
        cwtY = cwtX;
    else
        cwtY = getCWT(x(:, 2), x_len_adjusted, wft);
    end
    
    % ... and then multiply it with the corresponding conjugate
    Ps = cwtX .* conj(cwtY);

    % Optional smoothing. This is recommended for removing noise which 
    % appears after inverse Fourier transform. If smoothing is not 
    % performed, singular spectral matrices may appear in the function 
    % 'getCanonicalCoherenceW()' in '../coherence/getCanonicalCoherenceW.m'
    % due to nearly the same wavelet spectra
    Ps = smoothSpectrum(Ps, scales);
end

%% Auxiliaries
% Computation of the continuous wavelet transform of the data
function cwt = getCWT(x, x_len_adjusted, wft)
    x_len = size(x, 1);

    % We adjust the signal
    auxIndLeft = [x_len_adjusted : -1 : 1];
    auxIndRight = [x_len : -1 : x_len - x_len_adjusted + 1];
    xAdjusted = [conj(x(auxIndLeft)); x; conj(x(auxIndRight - 1))];

    % We compute the continuous Fourier transform of the adjusted time 
    % series ...
    xAdjustedFFT = fft(xAdjusted);
    % ... and compute the continuous wavelet transform of the adjusted time 
    % series via inverse Fourier transform
    cwt = ifft(repmat(transpose(xAdjustedFFT), size(wft, 1), 1) .* wft, [], 2);

    % Finally we keep only those values that correspond to the original 
    % time series' frequencies
    cwt = cwt(:, x_len_adjusted + 1 : x_len_adjusted + x_len);
end

% Computing the CWT of Morlet wavelet
function [wft, freqS, scales, coi] = MorletCWT(x_len, fs, preFreqs, x_adjustment, waveletSigma)
    % The list of dimensions of the quantities used:
    % [x_len] = [samples] = samples
    % [waveletSigma] = [omegaCutoff] = [freqT] = radians per second (angular frequencies)
    % [minScale] = [maxScale] = [scales] = relative units (dimensionless quantities)
    % [waveletFreq] = [preFreqs] = [minFreq] = [freqS] = [coi] = hertz, or cycles per second
    % [fs] = samples per second

    fpo = 12;                                               % Number of scale-related frequencies per octave. The greater this parameter, the more accurate the CWT
    a0 = 2.0 ^ (1.0 / fpo);                                 % Scale ratio (or geometric increment) between two successive voices
    
    kappaSigma = exp(-0.5 * waveletSigma ^ 2.0);
    cSigma = 1.0 / sqrt(1.0 + exp(-waveletSigma ^ 2.0) - 2.0 * exp(-0.75 * waveletSigma ^ 2.0));
    amplitude = cSigma / pi ^ 0.25;
    
    % We determine the cutoff angular frequency, in radians per second
    cutoffLevel = 0.1;
    waveletFourierTransform = @(omega)cutoffLevel - amplitude * (exp(-0.5 * (omega - waveletSigma) ^ 2.0) - kappaSigma * exp(-0.5 * omega ^ 2.0));  % Fourier transform of the Morlet wavelet minus the cutoff value
    
    omegaMax = sqrt(2.0 * 750.0) + waveletSigma;
    if (waveletFourierTransform(waveletSigma) > 0.0)
        omegaCutoff = omegaMax;
    else    
        omegaCutoff = fzero(waveletFourierTransform, [waveletSigma omegaMax]);
    end
    
    % We make a correction for the preset scale-related frequency (if
    % relevant)
    minScale = omegaCutoff / pi / 1.0;                      % Preliminary minimum admissible scale, in relative units (radians per second / radians per cycle / cycles per second)
    maxScale = x_len / (2.0 * sqrt(2.0));                   % Preliminary maximum admissible scale, in relative units (samples / samples)

    if (maxScale < minScale * a0)
        maxScale = minScale * a0;
    end

    waveletFreq = waveletSigma / (2.0 * pi);                % The frequency of the sine harmonic of the Morlet wavelet, a constant, in hertz (or in cycles per second)    
    minFreq = waveletFreq * 1.0 / maxScale * fs;            % Minimum admissible scale-related frequency, in hertz

    % So, we make the correction
    if (preFreqs(1) < minFreq)
        preFreqs(1) = minFreq;                              % The corrected preset scale-related frequency, in hertz
    end

    % Now we prepare the time-related frequencies (usually denoted in
    % scientific literature as 'omega'), in radians per second ...
    freqT = makeFourierTransformGrid(x_len, x_adjustment);
    
    nFreq = length(freqT);
    
    % ... and scales, in relative units
    preFreqs = 2.0 * pi * preFreqs / fs;                    % We convert the scale-related frequencies from hertz to radians per sample
    minScale = waveletSigma / preFreqs(2) / 1.0;            % The minimum and ...
    maxScale = waveletSigma / preFreqs(1) / 1.0;            % ... maximum admissible scales are updated but keep the same dimensions: relative units

    maxNumOctaves = log2(maxScale / minScale);
    scales = minScale * a0 .^ [0 : floor(maxNumOctaves * fpo)]';        % Scales, in relative units

    nScales = length(scales);

    % Finally, we compute the CWT of the Morlet wavelet on the determined 
    % time-related frequencies and scales
    wft = zeros(nScales, nFreq);
    for j = 1 : nScales
        eexp1 = -0.5 * (scales(j) * freqT - waveletSigma) .^ 2.0;
        eexp1 = eexp1 .* (freqT > 0.0);
        
        eexp2 = -0.5 * (scales(j) * freqT) .^ 2.0;
        eexp2 = eexp2 .* (freqT > 0.0);
        
        wft(j, :) = amplitude * sqrt(scales(j)) * (exp(eexp1) - kappaSigma * exp(eexp2)) .* (freqT > 0.0);
    end
    
    % As the last action, we compute the actual scale-related frequencies
    % ...
    freqS = waveletFreq * 1.0 ./ scales * fs;               % Actual scale-related frequencies, in hertz

    % ... and the cone of influence
    xl = ceil(x_len / 2.0);
    if (mod(x_len, 2) == 1)
        samples = [[1 : xl] [xl - 1 : -1 : 1]];
    else
        samples = [[1 : xl] [xl : -1 : 1]];
    end
    
    coi = waveletFreq * sqrt(2.0) ./ samples * fs * 1.0;    % In hertz
    
    % We truncate the excessive values of the cone of influence
    coi(coi > freqS(1)) = freqS(1);
    coi(coi < freqS(end)) = freqS(end);
end

% Adjusting signal length
function x_len_adjusted = adjustSignalLength(x_len)
    if (x_len <= 100000)
        x_len_adjusted = floor(x_len / 2.0);
    else
        x_len_adjusted = ceil(log2(x_len));
    end
end

% Making a time-related frequency grid for performing the Fourier
% transforms
function freqT = makeFourierTransformGrid(x_len, x_adjustment)
    % This is the size of the time-related frequency range
    N = x_len + 2 * x_adjustment;

    % This is the first half of the frequency range ...
    freqT = [1 : floor(N / 2.0)];
    freqT = freqT * (2.0 * pi) / N;
    % ... while these are the indices for covering the second half
    aux_inds = [floor((N - 1) / 2.0) : -1 : 1];

    % So, we make the frequency range, in radians per second
    freqT = [0.0 freqT -freqT(aux_inds)];
end

% Smoothing the spectrum
function mtrx = smoothSpectrum(mtrx, scales)
    % We are to filter along the time axis ...
    x_len = size(mtrx, 2);

    nfft = 2 ^ nextpow2(x_len);                     % We enlarge the number of samples to the closest power of two

    freqT = 2.0 * pi / nfft * [1 : floor(nfft / 2.0)];
    aux_inds = [floor((nfft - 1) / 2.0) : -1 : 1];
    freqT = [0.0 freqT -freqT(aux_inds)];           % Time-related frequencies, in radians per second

    for k = 1 : size(mtrx, 1)
        aux = exp(-0.5 * (scales(k) * freqT) .^ 2.0);
        smth = ifft(aux .* fft(mtrx(k, :), nfft));
        mtrx(k, :) = smth(1 : x_len);
    end
    
    % ... and then across the scales
    scalesToSmooth = 12;                            % The greater this parameter, the less accurate the spectrum
    flt = 1.0 / scalesToSmooth * ones(scalesToSmooth, 1);
    mtrx = conv2(mtrx, flt, 'same');
end
