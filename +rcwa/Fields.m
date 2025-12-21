function flds = Fields(x, y, z, coeffs, Layers, xj, sys)
    % this function calculates field componenets at specified spatial coordinates)
    % x,y,z: position to evaluate fields, same unit as lambda0
    % only 2 of 3 can be vectors(creates 2D cross section)

    % coeffs: from SolveStack (for a single polarization)
    % Layers: cell array of layers
    % xj: interface positions
    % sys: system structure from SystemSetup

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell
    % check inputs for cross sections
    % currently can not handle 3dplots efficently, must iterate
    if sum(cellfun(@length, {x, y, z}) > ones(1, 3)) > 2
        error('RCWA:Fields', 'only 2 of x,y,z can be vectors')
    end

    % function F=Fields(x,z,t,r,FA,Layer.kz,xj,sys)
    %%%Plots for 2D RCWA
    N = length(Layers) - 2;

    E0x = zeros(length(sys.kx), length(z));
    E0y = zeros(length(sys.kx), length(z));
    D0z = zeros(length(sys.kx), length(z));
    H0x = zeros(length(sys.kx), length(z));
    H0y = zeros(length(sys.kx), length(z));
    H0z = zeros(length(sys.kx), length(z));

    %% Top
    ind = find(z <= xj(1));

    if ~isempty(ind)
        phaseP = exp(1i * bsxfun(@times, Layers{1}.kz, (z(ind) - xj(1))));
        phaseM = exp(1i * bsxfun(@times, Layers{1}.kz, (xj(1) - z(ind))));

        E0 = Layers{1}.F{1} * bsxfun(@times, coeffs.fwd(:, 1), phaseP) + Layers{1}.F{1} * bsxfun(@times, coeffs.bkwd(:, 1), phaseM);
        H0 = Layers{1}.F{2} * bsxfun(@times, coeffs.fwd(:, 1), phaseP) - Layers{1}.F{2} * bsxfun(@times, coeffs.bkwd(:, 1), phaseM);

        E0x(:, ind) = E0(1:end / 2, :);
        E0y(:, ind) = E0(end / 2 + 1:end, :);
        H0x(:, ind) = H0(1:end / 2, :);
        H0y(:, ind) = H0(end / 2 + 1:end, :);

        D0z(:, ind) = (bsxfun(@times, sys.kx, H0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, H0(1:end / 2, :)))) / sys.omg;
        H0z(:, ind) = (bsxfun(@times, sys.kx, E0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, E0(1:end / 2, :)))) / sys.omg;

    end

    %% Middle
    if N ~= 0

        for nn = 2:N + 1
            ind = find(z > xj(nn - 1) & z <= xj(nn));

            if ~isempty(ind)
                phaseP = exp(1i * bsxfun(@times, Layers{nn}.kz, (z(ind) - xj(nn - 1))));
                phaseM = exp(1i * bsxfun(@times, Layers{nn}.kz, (xj(nn) - z(ind))));

                E0 = Layers{nn}.F{1} * bsxfun(@times, coeffs.fwd(:, nn), phaseP) + Layers{nn}.F{1} * bsxfun(@times, coeffs.bkwd(:, nn), phaseM);
                H0 = Layers{nn}.F{2} * bsxfun(@times, coeffs.fwd(:, nn), phaseP) - Layers{nn}.F{2} * bsxfun(@times, coeffs.bkwd(:, nn), phaseM);

                E0x(:, ind) = E0(1:end / 2, :);
                E0y(:, ind) = E0(end / 2 + 1:end, :);
                H0x(:, ind) = H0(1:end / 2, :);
                H0y(:, ind) = H0(end / 2 + 1:end, :);

                D0z(:, ind) = (bsxfun(@times, sys.kx, H0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, H0(1:end / 2, :)))) / sys.omg;
                H0z(:, ind) = (bsxfun(@times, sys.kx, E0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, E0(1:end / 2, :)))) / sys.omg;
            end

        end

    end

    %% End
    ind = find(z > xj(end));

    if ~isempty(ind)
        phaseP = exp(1i * bsxfun(@times, Layers{end}.kz, (z(ind) - xj(end))));

        E0 = Layers{end}.F{1} * bsxfun(@times, coeffs.fwd(:, end), phaseP);
        H0 = Layers{end}.F{2} * bsxfun(@times, coeffs.fwd(:, end), phaseP);

        E0x(:, ind) = E0(1:end / 2, :);
        E0y(:, ind) = E0(end / 2 + 1:end, :);
        H0x(:, ind) = H0(1:end / 2, :);
        H0y(:, ind) = H0(end / 2 + 1:end, :);

        D0z(:, ind) = (bsxfun(@times, sys.kx, H0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, H0(1:end / 2, :)))) / sys.omg;
        H0z(:, ind) = (bsxfun(@times, sys.kx, E0(end / 2 + 1:end, :) - bsxfun(@times, sys.ky, E0(1:end / 2, :)))) / sys.omg;
    end

    phiX = exp(1i * bsxfun(@times, [sys.kx], x));
    phiY = exp(1i * bsxfun(@times, [sys.ky], y));

    if length(x) == 1 || length(y) == 1
        phiXY = bsxfun(@times, phiX.', phiY.');

        Ex = phiXY * E0x;
        Ey = phiXY * E0y;
        Dz = phiXY * D0z;
        Hx = phiXY * H0x;
        Hy = phiXY * H0y;
        Hz = phiXY * H0z;
        flds = struct('Ex', Ex.', 'Ey', Ey.', 'Dz', Dz.', 'Hx', Hx.', 'Hy', Hy.', 'Hz', Hz.');
    else
        Ex = phiY.' * bsxfun(@times, E0x, phiX);
        Ey = phiY.' * bsxfun(@times, E0y, phiX);
        Dz = phiY.' * bsxfun(@times, D0z, phiX);
        Hx = phiY.' * bsxfun(@times, H0x, phiX);
        Hy = phiY.' * bsxfun(@times, H0y, phiX);
        Hz = phiY.' * bsxfun(@times, H0z, phiX);
        flds = struct('Ex', Ex, 'Ey', Ey, 'Dz', Dz, 'Hx', Hx, 'Hy', Hy, 'Hz', Hz);
    end

    % flds:
    % Ex, Ey, Dz: Electric field components(complex, with arbitrary units)
    % Hx, Hy, Hz : Magenetic field components(complex, with arbitrary units)
    % the choice of Hz and Dz reflects the most natural computationally efficient representation of Maxwell's equation in periodic, inhomogeneous media. media

end
