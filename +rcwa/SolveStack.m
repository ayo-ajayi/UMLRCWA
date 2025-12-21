function [Ref, Trans, coeffs] = SolveStack(layers, xj, a0, sys)
    % cell array of layer structures {layer1,layer2,...}
    %xj: z position of the interfaces between 0 and 1 interval
    % unit same as lambda0: example:
    % xj = [0, 0.1, 0.4, 0.7, 1.0]
    %       ↑    ↑    ↑    ↑    ↑
    %       │    │    │    │    └─ Interface 5: between Layer 5 & 6
    %       │    │    │    └────── Interface 4: between Layer 4 & 5
    %       │    │    └─────────── Interface 3: between Layer 3 & 4
    %       │    └──────────────── Interface 2: between Layer 2 & 3
    %       └───────────────────── Interface 1: between Layer 1 & 2

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell
    if ~iscell(layers)
        error('Something went very wrong')
    end

    Tmat = cell(length(layers) - 1, 1);

    for ii = 1:(length(layers) - 1)
        Tmat{ii} = TMatrix(layers{ii}.F, layers{ii + 1}.F);
    end

    kzArr = zeros(2 * sys.nTot, length(layers) - 2);

    for ii = 1:length(layers) - 2
        kzArr(:, ii) = layers{ii + 1}.kz;
    end

    [t, r, Cpl, Cmin] = ET(Tmat, kzArr, diff(xj), a0);
    [Ref, Trans] = DE(layers{1}.F, layers{end}.F, a0, r, t, 2);

    fwd = permute(cat(3, a0, Cpl, t), [1, 3, 2]);
    fwd = mat2cell(fwd, size(fwd, 1), size(fwd, 2), ones(size(a0, 2), 1));
    fwd = reshape(fwd, [size(fwd, 3), 1]);

    bkwd = permute(cat(3, r, Cmin, r * 0), [1, 3, 2]);
    bkwd = mat2cell(bkwd, size(bkwd, 1), size(bkwd, 2), ones(size(a0, 2), 1));
    bkwd = reshape(bkwd, [size(bkwd, 3), 1]);

    coeffs = struct('fwd', fwd, 'bkwd', bkwd);
    % coeffs=reshape(coeffs,[size(a0,2),1])
    % Ref: Reflection Coefficients (nTot x Npol) - dimensionless intensity for each diffraction order
    % Trans: Transmission Coefficients (nTot x Npol) - dimensionless intensity for each diffraction order
    % coeffs: Structure containing forward/backward propagating field coefficients at each layer
    % coeffs.fwd - travelling in +z direction
    % coeffs.bkwd - travelling in -z direction
    % energy conservation = sum(Ref)+sum(Trans) ≈ 1

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T] = TMatrix(F0, F1)
    %TMatrixMetaGrating
    %Finds the Tmatrix, for the following Matrix Problem
    %Solves the Following System of Equations
    % |W0  W0||i+| |W1  W1||t |
    % |      ||  |=|      ||  |
    % |V0 -V0||r | |V1 -V1||i-|
    %
    % Solution:
    % |i+| |T1 T2||t |
    % |  |=|     ||  |
    % |r | |T2 T1||i-|

    a = F0{1} \ F1{1};
    b = F0{2} \ F1{2};
    T{1} = (a + b) / 2;
    T{2} = (a - b) / 2;
end

function [t, r, Cpl, Cmin] = ET(T, kz, d, i0)
    %%%
    %caculate T & R through a stack, without any growing exponentials, or
    %inversion of decaying expoentials.
    %Must have at least one intermediate layer???
    %
    % Assume T is passed as a cellular array
    %                              |T{i}(1) T{i}(2)|
    % If T{i} is 2 elements,  T{i}=|               |
    %                              |T{i}(2) T{i}(1)|
    %
    %
    %                              |T{i}(1,1) T{i}(1,2)|
    % If T{i} is 4 elements,  T{i}=|                   |
    %                              |T{i}(2,1) T{i}(2,2)|

    % kz -> kz for each intermediate layer, [kz1,...,kzN]
    % d  -> height of each intermediate layer [d1,...dN]

    %i0 is a vector which corresponds to an incident wave
    %i0 may be a matrix of the form [i0(1),i0(2),...i0(N)];

    %Output: Cpl/Cmin will be off the shape (:,i0,N)
    %t,r will be of the shape(:,i0)

    %This version will calculate all intermediate coefficients according to the
    %following conventions

    %phase                exp(i*kz*z)-->
    %                |        |   |          |
    %            i0->|c+_1->  |   |c+_N-1->  |t->
    %                |        |...|          |
    %                |        |   |          |
    %             <-r|  <-c-_1|   |  <-c-_N-1|
    %                |        |   |          |
    %interface#      1        2  N-1         N
    %layer #      0      1             N         N+1

    %find total number of interfaces
    N = length(T);
    Ni = size(i0, 2);

    if isempty(kz)
        nKz = length(i0(1, :));
        kz = 0;
        d = 0;
    else
        nKz = length(kz(:, 1));
    end

    %initialize storage matricies/cells
    Cpl = zeros(nKz, Ni, N - 1);
    Cmin = Cpl;
    a = cell(N, 1);
    b = a;

    %If there is only a single interface then we are essentially done
    if N == 1

        if length(T{1}) == 2
            t = T{1}{1} \ i0;
            r = T{1}{2} * t;
        else
            t = T{1}{1, 1} \ i0;
            r = T{1}{2, 1} * t;
        end

        return
    end

    phi = exp(bsxfun(@times, kz, 1i * d));

    %Last Interface
    if length(T{1}) == 2
        a{N} = T{N}{1};
        b{N} = T{N}{2};
    else
        a{N} = T{N}{1, 1};
        b{N} = T{N}{2, 1};
    end

    Sn = (phi(:, N - 1) * phi(:, N - 1).') .* (b{N} / a{N});
    %iterate backwards throught the interfaces
    for n = (N - 1):-1:2

        if length(T{n}) == 2
            a{n} = T{n}{1} + T{n}{2} * Sn;
            b{n} = T{n}{2} + T{n}{1} * Sn;
        else
            a{n} = T{n}{1, 1} + T{n}{1, 2} * Sn;
            b{n} = T{n}{2, 1} + T{n}{2, 2} * Sn;
        end

        Sn = (phi(:, n - 1) * phi(:, n - 1).') .* (b{n} / a{n});
    end

    %Initial Interface
    if length(T{1}) == 2
        a{1} = T{1}{1} + T{1}{2} * Sn;
        b{1} = T{1}{2} + T{1}{1} * Sn;
    else
        a{1} = T{1}{1, 1} + T{1}{1, 2} * Sn;
        b{1} = T{1}{2, 1} + T{1}{2, 2} * Sn;
    end

    %%%Arg 1d exception

    %calculate all the intermediate fields
    Cpl(:, :, 1) = a{1} \ i0;
    r = b{1} * Cpl(:, :, 1);

    %loop back again through intermediate layers to find all coefficients
    for n = 2:N - 1
        Cpl(:, :, n) = a{n} \ bsxfun(@times, phi(:, n - 1), Cpl(:, :, n - 1));
        Cmin(:, :, n - 1) = b{n} * Cpl(:, :, n);
    end

    %final layer
    t = a{N} \ bsxfun(@times, phi(:, N - 1), Cpl(:, :, N - 1));
    Cmin(:, :, N - 1) = b{N} * t;
end

function [R, T, Sinc] = DE(F0, FN, a0, Cmin, Cpl, dim)
    %Calculates Diffraction Efficiency
    %dim (optional) -> 1,2 for (1D/2D) (defaults to 1D)

    if nargin == 5
        dim = 1;
    end

    if dim == 1
        Sinc = sum(real((F0{1} * a0) .* conj(F0{2} * a0)));
        R = real((F0{1} * Cmin) .* conj(F0{2} * Cmin) / Sinc);
        T = real((FN{1} * Cpl) .* conj(FN{2} * Cpl) / Sinc);
    elseif dim == 2
        Sinc = sum(real((F0{1}(1:end / 2, :) * a0) .* conj(F0{2}((end / 2 + 1):end, :) * a0))- ...
            real((F0{1}((end / 2 + 1):end, :) * a0) .* conj(F0{2}(1:end / 2, :) * a0)), 1);

        R = real(bsxfun(@times, (F0{1}(1:end / 2, :) * Cmin) .* conj(F0{2}(end / 2 + 1:end, :) * Cmin)- ...
            (F0{1}((end / 2 + 1):end, :) * Cmin) .* conj(F0{2}(1:end / 2, :) * Cmin), 1 ./ Sinc));

        T = real(bsxfun(@times, (FN{1}(1:end / 2, :) * Cpl) .* conj(FN{2}((end / 2 + 1):end, :) * Cpl)- ...
            (FN{1}((end / 2 + 1):end, :) * Cpl) .* conj(FN{2}(1:end / 2, :) * Cpl), 1 ./ Sinc));
    end

end
