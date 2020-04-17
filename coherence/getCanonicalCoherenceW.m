function [evt, ev, freq] = getCanonicalCoherenceW(xx)
% [evt, ev, freq] = getCanonicalCoherenceW(xx)
%
% This code implements the canonical coherence analysis of multivariate 
% data. The computation is performed using a wavelet spectrum estimation. 
% For details, please refer to [1].
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
%
% The end user is granted perpetual permission to reproduce, adapt, and/or 
% distribute this code, provided that an appropriate link is given to the 
% original repository it was downloaded from.

    %% Auxiliaries
    %tic;    % DEBUG
    
    N = size(xx, 2);    % Number of variates in the vector time series

    % Determining the frequency range
    [~, aux, freq] = wcoherence(xx(:, 1), xx(:, 1), 1.0);
    freq_len = length(freq);

    t_len = size(aux, 2);

    %% Computing
    ev = zeros(freq_len, t_len, N);
    for ic = 1 : N
        %% We split the original N-variate signal into an (N - 1)-variate 
        % and a single-variate ones
        x = xx(:, [1 : ic - 1 ic + 1 : N]);     % (N - 1)-variate (i.e. vector) time series
        y = xx(:, ic);                          % Single-variate (i.e. scalar) time series

        %% We compute the spectral matrices ...
        % ... for the first signal, ...
        Sxx = zeros(N - 1, N - 1, freq_len, t_len);
        for i = 1 : N - 1
            for j = 1 : N - 1
                [~, Sxx(i, j, :, :), ~] = wcoherence(x(:, i), x(:, j), 1.0);
            end
        end

        % ... for the mixture of the first and second signals, ...
        Sxy = zeros(N - 1, 1, freq_len, t_len);
        for i = 1 : N - 1
            [~, Sxy(i, 1, :, :), ~] = wcoherence(x(:, i), y(:, 1), 1.0);
        end

        Syx = zeros(1, N - 1, freq_len, t_len);
        for j = 1 : N - 1
            [~, Syx(1, j, :, :), ~] = wcoherence(y(:, 1), x(:, j), 1.0);
        end

        % ... and for the second signal
        Syy = zeros(1, 1, freq_len, t_len);
        [~, Syy(1, 1, :, :), ~] = wcoherence(y(:, 1), y(:, 1), 1.0);

        %% We compute the matrix whose eigenvalues are PCA-related measures 
        % of coherence and keep the maximum value
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

        aux1 = cellfun(@mtimes, Sxx_aux, Sxy_aux, 'UniformOutput', false);
        aux2 = cellfun(@mtimes, aux1, Syy_aux, 'UniformOutput', false);
        UU = cellfun(@mtimes, aux2, Syx_aux, 'UniformOutput', false);
        
        aux = cellfun(@eig, UU, 'UniformOutput', false);
        aux = cellfun(@real, aux, 'UniformOutput', false);
        ev(:, :, ic) = reshape(cellfun(@max, aux), freq_len, t_len);
        %}

        % Normalisation
        ev(:, :, ic) = ev(:, :, ic) / max(max(abs(ev(:, :, ic))));
    end

    % Finally we compute the total coherence frequency-wise
    evt = sqrt(prod(ev, 3));

    %toc;    % DEBUG
end

% This auxiliary function reshapes a 4-D array into a sequence of 2-D
% matrices
function res = toMatrixArray(mtrx)
    aux = reshape(mtrx, size(mtrx, 1), size(mtrx, 4) * size(mtrx, 3) * size(mtrx, 2)); 
    res = mat2cell(aux, size(mtrx, 1), ones(1, size(mtrx, 4) * size(mtrx, 3)) * size(mtrx, 2));
end    
