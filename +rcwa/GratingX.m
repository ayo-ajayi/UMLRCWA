% Copyright (c) 2014-2017 by University of Massachusetts Lowell
function [layer] = GratingX(epsOut, epsIn, len, sys)
    % epsOut: permittivity outside the grating lines(dimensionless)
    % epsIn: permittivity inside the grating lines(dimensionless)
    % len: width of grating lines(same unit as Lam from SystemSetup)
    %sys: system structure from SystemSetup
    % One Period (Lam)
    %    |<--------------->|

    %    ┌────┬────────────┬────┬────────────┬────┬──
    %    │    │            │    │            │    │
    %    │ In │    Out     │ In │    Out     │ In │  ... (repeats)
    %    │    │            │    │            │    │
    %    └────┴────────────┴────┴────────────┴────┴──
    %    |<-->|<--------->|
    %     len   Lam - len

    [epsM, epsInvM] = GratingX2D([epsOut, epsIn], len, sys.Lam, sys.maxM);
    NV = NormalVectorField_GratingX2D(sys.maxM);
    [F, kz] = Dispersion2D(epsM, epsInvM, sys, NV);
    layer = struct('F', {F}, 'kz', kz);
    %same return from homogeneous layer but for grating
    % len/Lam(1) fill factor from o to 1
end

function [epsMat, epsInvMat] = GratingX2D(epsi, len1, Lam, maxM)
    %Given the indices and geometry out puts eps toeplitz matrix

    % Create eps convolution index Matrix
    NM = 2 * maxM + 1;
    nTot = prod(NM);
    %xInd=(mod((1:nTot)-1,NM(1))-(NM(1)-1)/2);
    %yInd=(ceil(((1:nTot))/NM(1))-(NM(2)+1)/2);

    xInd = mod((1:nTot) - 1, NM(1)) + 1;
    yInd = ceil((1:nTot) / NM(1));

    xInd = bsxfun(@minus, xInd', xInd);
    yInd = bsxfun(@minus, yInd', yInd);

    q0 = 2 * pi / Lam(1);

    if length(epsi) == 1
        eps1 = epsi;
        eps2 = eps1;
    else
        eps1 = epsi(1);
        eps2 = epsi(2);
    end

    %Construct normal epsilon matrix
    tempX = zeros(size(xInd));
    tempY = zeros(size(yInd));

    tempX = 2 * (eps1 - eps2) * sin(len1 * q0 * xInd / 2) ./ (q0 * xInd) / Lam(1);
    tempX(xInd == 0) = (eps1 * len1 + eps2 * (Lam(1) - len1)) / Lam(1);
    tempY(yInd == 0) = 1;
    epsMat = tempX .* tempY;

    tempX = zeros(size(xInd));
    tempY = zeros(size(yInd));

    tempX = 2 * (1 / eps1 - 1 / eps2) * sin(len1 * q0 * xInd / 2) ./ (q0 * xInd) / Lam(1);
    tempX(xInd == 0) = (1 / eps1 * len1 + 1 / eps2 * (Lam(1) - len1)) / Lam(1);
    tempY(yInd == 0) = 1;
    epsInvMat = tempX .* tempY;
end

function NV = NormalVectorField_GratingX2D(maxM)
    O = diag(ones(1, prod(2 * maxM + 1)));
    Z = zeros(prod(2 * maxM + 1));
    NV = {Z, Z; Z, O};
end

function [F, kzP] = Dispersion2D(epsM, epsInvM, sys, NV)
    %This constructs and solves the RCWA eigenvalue problem
    %Returns: kzArr eigenvectors Forward Moving Waves Only!
    %smallNum=1e-10;

    if nargin == 4%Normal Vector Field
        NV{1} = epsM * NV{1} - epsInvM \ NV{1}; %xx
        NV{2} = epsM * NV{2} - epsInvM \ NV{2}; %yx
        NV{3} = epsM * NV{3} - epsInvM \ NV{3}; %xy
        NV{4} = epsM * NV{4} - epsInvM \ NV{4}; %yy
    elseif nargin == 2%Moriham
        sys = epsInvM;
        epsInvM = epsM;
        NV = {0, 0; 0, 0};
    elseif nargin == 3%Li
        NV = {0, 0; 0, 0};
        NV{1} = epsM - epsInvM \ eye(size(epsInvM));
    end

    %Create idenity Matrix for easy calc
    Idty = eye(size(epsM));

    %Normalize
    nx = sys.kx / sys.omg;
    ny = sys.ky / sys.omg;

    % Now we can construct the eigenvalue problem
    P = [(nx * ny.') .* (epsM \ Idty), Idty - (nx * nx.') .* (epsM \ Idty); ...
            (ny * ny.') .* (epsM \ Idty) - Idty, -(ny * nx.') .* (epsM \ Idty)];

    Q = [-diag(nx .* ny) + NV{3}, diag(nx.^2) - epsM + NV{4}; ...
            epsM - diag(ny.^2) - NV{1}, diag(ny .* nx) - NV{2}];

    F = cell(2, 1);
    %Now we reduce the system by factor of 1/2 by solving the wave equation
    %Since eigen value analysis scales as N^3 this is a net gain of x8.
    [F{1}, kzP] = eig(P * Q);

    %Eigenvectors are not necessarily normalized due to 'nobalance' if needed
    %[F{1},kzP]=eig(P*Q,'nobalance');
    %F{1}=bsxfun(@mrdivide,F{1},sqrt(sum(abs(F{1}).^2)));

    kzP = diag(kzP);

    %Now we need to take the square root of the eigenvalues
    %And form the magnetic eigenvecots
    %kzP=mySqrt(kzP,-pi/2);
    kzP = sqrt(kzP);

    %%Important clever branch cut do not prevent the sign on imag(kzP) from
    %%being all the same, so since we know the values are +/-, and we want the
    %%fields to decay the same, we flip the propegation, and hence the H-field
    %%will flip below!
    kzP(imag(kzP) < 0)=-kzP(imag(kzP) < 0);

    F{2} = bsxfun(@times, Q * F{1}, (1 ./ kzP).');
    kzP = kzP * sys.omg;

end
