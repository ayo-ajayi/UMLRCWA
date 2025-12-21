function [layer] = HomogenousLayer(epsilon, sys)
    % creates eigenmode solutions for a homogeneous (uniform) layer
    % epsilon: permittivity(dimensionless): scalar(isotropic material) or [epsXY, epsZ] (uniaxial anisotropic material)
    %sys: syste structure from SystemSetup

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    %% Calculates eigenFields and kz+ for a uniaxial materials

    if length(epsilon) == 1
        epsXY = epsilon;
        epsZ = epsilon;
    elseif length(epsilon) == 2
        epsXY = epsilon(1);
        epsZ = epsilon(2);
    else
        error(['epsilon must be a 1 or 2 component array.  epsilon=[epsIsotropic] or epsilon=[epsXY,epsZ]'])
    end

    nx = sys.kx / sys.omg;
    ny = sys.ky / sys.omg;

    nzTE = sqrt(epsXY - (nx.^2 + ny.^2));
    nzTM = sqrt(epsXY - (nx.^2 + ny.^2) * epsXY / epsZ);

    norm = sqrt((nx.^2 + ny.^2));

    F = cell(2, 1);

    F{1} = [diag(-ny ./ norm), diag(nzTM .* nx / epsXY ./ norm); ...
            diag(nx ./ norm), diag(nzTM .* ny / epsXY ./ norm)];

    F{2} = [diag(-nx .* nzTE ./ norm), diag(-ny ./ norm); ...
            diag(-ny .* nzTE ./ norm), diag(nx ./ norm)];

    %special case for normal incidence.
    ind = (nx == 0 & ny == 0);

    if sum(ind) == 1
        F{1}([ind; ind], [ind; ind]) = [sind(sys.phi), cosd(sys.phi) * nzTM(ind) / epsXY; ...
                                    cosd(sys.phi), -sind(sys.phi) * nzTM(ind) / epsXY];
        F{2}([ind; ind], [ind; ind]) = [-nzTE(ind) * cosd(sys.phi), sind(sys.phi); ...
                                    nzTE(ind) * sind(sys.phi), cosd(sys.phi)];

    end

    kz = sys.omg * [nzTE; nzTM];

    layer = struct('F', {F}, 'kz', kz);
    % layer: F, kz
    % F:[F1=electric field eigenvectors(TM modes)-(2xnTot)x(2xnTot), F2=magnetic field eigenvectors(TE modes)-(2xnTot)x(2xnTot)]
    % kz: z-component wavevectors for all modes, size: (2xnTot)x1; first half: TE modes, second half: TM modes
end
