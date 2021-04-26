function [evt, ev, freq] = getCanonicalCoherence(ddx, fs, isInfo)
% [evt, ev, freq] = getCanonicalCoherence(ddx, fs, isInfo)
%
% This code implements canonical coherence analysis of multivariate 
% data. The computation is performed using a non-parametric Fourier 
% spectrum estimation. For details, please refer to [1, 2].
%
% At the input, 'ddx' is a multivariate time series of the second 
% increments of the physical observable 'ddx(t) = (ddx_1(t), ...,
% ddx_N(t))', 'fs' is the sample rate (optional), 'isInfo' is the flag
% prescribing to output a message for each variate processed (optional).
%
% At the output, 'evt' is the total coherence coefficient over the 
% frequency range, 'ev' is the array of the partial coherence coefficients 
% over the frequency range, 'freq' is the frequency range over which the 
% coherence spectrum has been computed.
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
    if (nargin == 1)
        fs = 1.0;
        isInfo = false;
    elseif (nargin == 2)
        isInfo = false;
    end
    
    N = size(ddx, 2);    % Number of variates in the vector time series

    % Determining the frequency range
    [~, freq] = cpsd(ddx(:, 1), ddx(:, 1), [], [], [], fs);
    freq_len = length(freq);

    %% Computing
    ev = zeros(freq_len, N);
    for ic = 1 : N
        %% We split the original N-variate signal into an (N - 1)-variate 
        % and a single-variate ones
        x = ddx(:, [1 : ic - 1 ic + 1 : N]);     % (N - 1)-variate (i.e. vector) time series
        y = ddx(:, ic);                          % Single-variate (i.e. scalar) time series

        %% We compute the spectral matrices (see [1, p. 111]) ...
        % ... for the first signal, ...
        Sxx = zeros(N - 1, N - 1, freq_len);
        for i = 1 : N - 1
            for j = 1 : N - 1
                Sxx(i, j, :) = cpsd(x(:, i), x(:, j), [], [], [], fs);
            end
        end

        % ... for the mixture of the first and second signals, ...
        Sxy = zeros(N - 1, 1, freq_len);
        for i = 1 : N - 1
            Sxy(i, 1, :) = cpsd(x(:, i), y(:, 1), [], [], [], fs);
        end

        Syx = zeros(1, N - 1, freq_len);
        for j = 1 : N - 1
            Syx(1, j, :) = cpsd(y(:, 1), x(:, j), [], [], [], fs);
        end

        % ... and for the second signal
        Syy = zeros(1, 1, freq_len);
        Syy(1, 1, :) = cpsd(y(:, 1), y(:, 1), [], [], [], fs);

        %% We compute the matrix, whose eigenvalues are PCA-related 
        % measures of coherence, and keep the maximum value
        %{
        % Way 1 (element-wise)
        for i = 1 : freq_len
            % We implement formula (2.3.1) from [1]
            U = inv(sqrtm(Sxx(:, :, i))) * Sxy(:, :, i) * inv(Syy(:, :, i)) * Syx(:, :, i) * inv(sqrtm(Sxx(:, :, i)));

            ev(i, ic) = max(real(eig(U)));
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
        ev(:, ic) = cellfun(@max, aux);
        %}
        
        % Informational message to the standard output
        if (isInfo)
            fprintf('CCA: variate %d of %d processed\n', ic, N);
        end
    end

    % Finally we compute the total coherence frequency-wise
    evt = prod(ev, 2) .^ (1.0 / N);
end

% This auxiliary function reshapes a 3-D array into a sequence of 2-D
% matrices
function res = toMatrixArray(mtrx)
    aux = reshape(mtrx, size(mtrx, 1), size(mtrx, 3) * size(mtrx, 2)); 
    res = mat2cell(aux, size(mtrx, 1), ones(1, size(mtrx, 3)) * size(mtrx, 2));
end    
