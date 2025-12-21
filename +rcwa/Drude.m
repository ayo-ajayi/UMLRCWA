function [epsilon] = Drude(omg, eps_b, lam_p, G)
    % Copyright (c) 2014-2017 by University of Massachusetts Lowell
    %DRUDE This function calculates permittivity from the
    % parameters given by UIUC, for doped dielectrics.
    % INPUT
    % omg: Angular Frequency [rad/um] (can be an array)
    % eps_b: bound permittiviy
    % lam_p: Plasma Wavelength [um]
    % G: Scattering Rate [Hz]
    % OUTPUT
    % epsilon: permitivity
    wp = 2 * pi / lam_p;
    %lam0=2*pi/omg;
    g = G / 3e14;
    epsilon = eps_b .* (1 - wp^2 ./ (omg.^2 + 1i * g * omg));
end
