function [a0, i0] = IncidentLight(psi, maxM)
    % creates the incident wave excitation
    % psi: polarization angle in degree
    % 0 = TE, 90 = TM
    % Can be vector for both polarizations
    % maxM: diffraction orders

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    %a0 vector excites (0,0)th order mode
    % excitation coefficient vector (2xnTotxNpol); where Npol is length(psi)
    %i0 returns location of (0,0)th mode

    if nargin == 1
        %psi will actually be maxM
        a0 = zeros(2 * psi + 1, 1);
        i0 = psi + 1;
        a0(i0) = 1;
    else
        nTot = prod(2 * maxM + 1);
        a0 = zeros(2 * nTot, length(psi));
        i0 = (nTot + 1) / 2;
        a0([0, nTot] + i0, :) = [cosd(psi); sind(psi)];
    end

end
