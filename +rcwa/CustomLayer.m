% Copyright (c) 2014-2017 by University of Massachusetts Lowell
function [layer] = CustomLayer(A, sys)
    % Allows any permittivity pattern
    % A is a permittivity_map: 2D matrix defining ε(x,y)
    % sys: system structure from SystemSetup
    %this function is very flexible but may converge slowly
    [epsM, epsInvM] = rcwa.convmat(A, sys.maxM);
    [F, kz] = Dispersion2D(epsM, epsInvM, sys);
    layer = struct('F', {F}, 'kz', kz, 'A', A);

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

    %% Now we can construct the eigenvalue problem
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
