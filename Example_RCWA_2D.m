clear
clc

import rcwa.*

%% System Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lam0 = 1;
omg = 2 * pi / lam0; %k=2pi/lambda

theta = 45;
phi = 0;

Lam = [1.75, 1.75]; % period X,Y
maxM = 7 .* [1, 1]; % number of orders
sys = SystemSetup([theta, phi], omg, Lam, maxM);
disp(sys);

%% Layer 1 - Superstrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Layer{1} = HomogenousLayer(1, sys);

%% Layer 2 - Grating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsOut = 1; % (Lam(1)-len)/Lam(1)
epsIn=-10 + 1i; % len/Lam(1)
len = .7;
Layer{2} = GratingX(epsOut, epsIn, len, sys);

%% Layer 3 - Disk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsIn = 1;
epsOut=-10 + 1i;
R = .5; %Radius
Nsp = 2^10 * [1, 1]; %Spatial grid used to approximate Disk%1024 x 1024 grid
Layer{3} = Disk(epsOut, epsIn, R, Nsp, sys);

% Plot the permittivity of the disk
% (Not implemented for gratings or homogenous layer)
% figure()
% pcolor(real(Layer{3}.A)),shading flat

%% Layer 4 - Square %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsIn = 1;
epsOut=-10 + 1i;
L = .5; %Length of a Side
Nsp = 2^10 * [1, 1]; %Spatial grid used to approximate Disk
%Layer{4} = Square(epsOut,epsIn,R,Nsp,sys); %original code
Layer{4} = Square(epsOut, epsIn, L, Nsp, sys); %original code

% Plot the permittivity of the disk
%(Not implemented for gratings or homogenous layer)
% figure()
% pcolor(real(Layer{4}.A)),shading flat

%% Layer 5 - Custom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsIn = 1;
epsOut=-10 + 1i;
Nsp = 2^10 * [1, 1];
R = .25;
N = 3;

% Create a squircle |x|^N+|y|^N=|R|^N
xv = linspace(-Lam(1) / 2, Lam(2) / 2, Nsp(1));
yv = linspace(-Lam(1) / 2, Lam(2) / 2, Nsp(2));
[xArr, yArr] = meshgrid(xv, yv);
Permittivity_Map = ones(size(xArr)) * epsOut;
ind = (abs(xArr).^N + abs(yArr).^N <= abs(R).^N);
Permittivity_Map(ind) = epsIn;

% figure()
% pcolor(real(Permittivity_Map)),shading flat

% Perform RCWA on custom permittivoty map
% may weakly convergence as function of Nsp,maxM
Layer{5} = CustomLayer(Permittivity_Map, sys);

%% Layer 6 - Substrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (same as first layer, no reason to recalculate)
Layer{6} = Layer{1};

%% Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xj = [0, .1, .4, .7, 1]; %location of interfaces [(# of layers)-1 entries]
psi = [0, 90]; %TE=0%TM=90
[a0] = IncidentLight(psi, maxM); % creates inital excitation on superstrate
[Ref, Trans, coeffs] = SolveStack(Layer, xj, a0, sys);

%% View Reflection/Transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table of selected diffraction orders
orders = [0, -1
    -1, 0;
    0, 0;
    1, 0
    0, 1];
ind = ord2ind(orders, sys);
table(orders, Ref(ind, 1), Ref(ind, 2), 'VariableNames', {'R_nm', 'TE', 'TM'})

% Table of all diffraction orders
ind = 1:sys.nTot;
orders = ind2ord(ind, sys);
table(orders, Ref(ind, 1), Trans(ind, 2), 'VariableNames', {'T_nm', 'TE', 'TM'})

%% Fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently calculates cross-sections only
x=-Lam(1) / 2:.01:Lam(1) / 2;
y = 0;
z = (xj(1) - lam0):.01:(xj(end) + lam0);

% Fields for TE, xz-cross section
TEFlds = Fields(x, y, z, coeffs(1), Layer, xj, sys);
figure()
subplot(2, 3, 1)
imagesc(x, z, abs(TEFlds.Ex))
subplot(2, 3, 2)
imagesc(x, z, abs(TEFlds.Ey))
subplot(2, 3, 3)
imagesc(x, z, abs(TEFlds.Dz))
subplot(2, 3, 4)
imagesc(x, z, abs(TEFlds.Hx))
subplot(2, 3, 5)
imagesc(x, z, abs(TEFlds.Hy))
subplot(2, 3, 6)
imagesc(x, z, abs(TEFlds.Hz))

% Fields for TM, xy-cross section
x=-Lam(1) / 2:.01:Lam(1) / 2;
y=-Lam(2) / 2:.05:Lam(2) / 2;
z = 0;
TMFlds = Fields(x, y, z, coeffs(2), Layer, xj, sys);
figure()
subplot(2, 3, 1)
imagesc(x, y, real(TMFlds.Ex))
subplot(2, 3, 2)
imagesc(x, y, real(TMFlds.Ey))
subplot(2, 3, 3)
imagesc(x, y, real(TMFlds.Dz))
subplot(2, 3, 4)
imagesc(x, y, real(TMFlds.Hx))
subplot(2, 3, 5)
imagesc(x, y, real(TMFlds.Hy))
subplot(2, 3, 6)
imagesc(x, y, real(TMFlds.Hz))

figs = findall(0, 'Type', 'figure');

for i = 1:length(figs)
    waitfor(figs(i));
end

clearvars;

%  matlab -batch "Example_RCWA_2D"
