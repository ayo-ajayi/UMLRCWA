function [sysStruct] = SystemSetup(angles, omg, Lam, maxM)
    %anges: [theta, phi] incident angles in degrees
    %theta: polar angle(0=normal incidence)
    %phi: azimuthal angle (0=plane of incidence is x-z plane)
    %omg: angular freq. 2pi/lambda0
    % Lam:[Lx,Ly] periods in x and y directions (same unit as lambda0)
    % maxM: [Mx, My] max diffraction orders in each direction (integer values)

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell
    %Angles-->[theta,phi]
    %Returns transverse components of the wavevector
    %If system is periodic supply Lam=[LamX,LamY], or Lam=[Lam,Lam]
    %For non-periodic system last two arguments are ignored.

    smallNum = 0; eps(1);

    if length(angles) == 1
        theta = angles(1);
        phi = 0;
    elseif length(angles) == 2
        theta = angles(1);
        phi = angles(2);
    end

    % if length(maxM)==1
    %     maxM=[maxM,maxM];
    % end

    if nargin == 2
        Lam = 1;
        maxM = [0, 0];
    end

    NM = 2 * maxM + 1;
    nTot = prod(NM);

    if length(Lam) == 1
        q0 = 2 * pi ./ [Lam, Lam];
    else
        q0 = 2 * pi ./ Lam;
    end

    kx0 = omg * sind(theta) * cosd(phi);
    kxm = kx0 + q0(1) * (mod((1:nTot) - 1, NM(1)) - maxM(1))';

    if length(maxM) == 2
        ky0 = omg * sind(theta) * sind(phi);
        kym = ky0 + q0(2) * (ceil(((1:nTot)) / NM(1)) - maxM(2) - 1)';
    else
        kym = 0;
    end

    sysStruct = struct('kx', kxm, 'ky', kym, 'theta', theta, 'phi', phi, 'omg', omg, 'Lam', Lam, 'maxM', maxM, 'smallNum', smallNum, 'nTot', nTot);
    %sysSTruct returns transverse components
    % kx,ky: transverse wavevector components for all diffraction orders (units: omgxrefractive index)
    % nTot: total number of diffraction orders=(2Mx+1)x(2My+1)

    % all variables are normalized to lambda0
    %lambda0 = 1 even tho lambda0 = 500nm for example
    % then Lam=[1.75, 1.75] means Lamx=1.75lambda0=875nm

    % Example:

    % lambda_real = 600e-9;  % meters
    % R_real = 0.3 * lambda_real = 180 nm
    % Lam_real = 1.5 * lambda_real = 900 nm

    % becomes

    % lam0 = 1;
    % R = 0.3;           % Normalized
    % Lam = [1.5, 1.5];  % Normalized

end
