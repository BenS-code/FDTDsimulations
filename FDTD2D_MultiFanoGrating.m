%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINITE DIFFERENCE TIME DOMAIN SIMULATION - 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
The 2D simulation (upgraded from Raymond Rampf's lectures) solves the
reflection & transmission spectrum of a pulse
impinges perpendicularly on a slab with material properties as summarized
in the pdf file in folder "FDTD_equations".

Spatial grid is based on "YEE grid" 

The results are for the TE and TM modes of the pulse 

Material properties, frequency range and the number of layers can be
modified.
%}

%% MULTI FANO RESONANCE GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZE MATLAB

close all; clear all; clc;

%% INITIALIZE FIGURE 

scrnsize = get(0,'screensize');
fig1 = figure('numbertitle','off','name','FDTD Simulation',...
    'color','white','position',scrnsize);

%% INITIALIZE COLORMAP
Cres = 256;
Gres = ceil(Cres/2);
Ccont = 2;
CMAP = zeros(Cres+Gres,3);
c1 = [0 0 1];
c2 = [1 1 1];
c3 = [1 0 0];
for nc = 1:ceil(Cres/2)
    f = (nc-1)/ceil(Cres/2);
    c = (1 - f^Ccont)*c1 + f^Ccont*c2;
    CMAP(nc,:) = c;
    c = (1 - f^(1/Ccont))*c2 + f^(1/Ccont)*c3;
    CMAP(ceil(Cres/2)+nc,:) = c;
    CMAP(ceil(Cres/2)+Gres+nc,:) = (1-f)*c2;
end

%% UNITS

m = 1;
cm = 1e-2 * m;
mm = 1e-3 * m;
um = 1e-6 * m;
nm = 1e-9 * m;

s = 1;
ms = 1e-3 * s;
us = 1e-6 * s;
ns = 1e-9 * s;
ps = 1e-12 * s;

Hz = 1;
kHz = 1e3 * Hz;
MHz = 1e6 * Hz;
GHz = 1e9 * Hz;
THz = 1e12 * Hz;
PHz = 1e15 * Hz;

F = 1;                  
H = 1;
Ohm = 1;

%% CONSTANTS

e0 = 8.8541878176 * 1e-12 * F/m;
u0 = 1.2566370614 * 1e-6 * H/m;
c0 = 1/sqrt(e0*u0) * m/s;
N0 = sqrt(u0/e0) * Ohm;

%% SIMULATION PROPERTIES

% OPTS.emax = 0.01;

% FREQUENCY RANGE
fmin = 0.3 * PHz;                    % minimum frequency
fmax = 0.5 * PHz;                    % maximum frequency
NFREQ = 300;                         % number of frequencies
FREQ = linspace(fmin,fmax,NFREQ);    % frequency array
        
% MATERIAL PROPERTIES
Wp_E1 = 0*13.72e15;
Gamma_E1 = 0*[0.08052 0.3661 0.5241].*1e15;
W0_E1 = 0*[0 0.6305 1.261].*1e15;
f_E1 = 0*[0.76 0.024 0.01];
er1_inf = 14;

Wp_M1 = 0;
Gamma_M1 = [0 0 0];
W0_M1 = [0 0 0];
f_M1 = [0 0 0];
ur1_inf = 1;

Wp_E2 = 0*13.72e15;
Gamma_E2 = 0*[0.08052 0.3661 0.5241].*1e15;
W0_E2 = 0*[0 0.6305 1.261].*1e15;
f_E2 = 0*[0.76 0.024 0.01];
er2_inf = 14;

Wp_M2 = 0;
Gamma_M2 = [0 0 0];
W0_M2 = [0 0 0];
f_M2 = [0 0 0];
ur2_inf = 1;

n_res = 1;%length(f_E1);

update_plot = 50;           % update plot every "update_plot" steps

% SLAB PROPERTIES
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  er_bg,ur_bg                                               background
%                        ________      ______     ____            _
%  er1,ur1              |        |    |      |   |    |
%                       |<--a1-->|    |<-b1->|   |<c1>|        layer1
%                       |________|____|______|___|____|_          _
%  er2,ur2              |             |          |      |
%                       |             |          |      |      layer2
%                       |_____________|__________|______|         _
%  er_sub,ur_sub        |             |          |      |
%                       |<----A1----->|<---B1--->|<-C1->|     substrate    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = 195 * nm;
a1 = 175 * nm;
B1 = 216 * nm;
b1 = 65 * nm;
C1 = 0 * nm;
c1 = 0 * nm;
layer1 = 20 * nm;
layer2 = 100 * nm;
er_bg = 1.0;
ur_bg = 1.0;
er_sub = 3.02;
ur_sub = 1.0;
er1 = er1_inf;
ur1 = ur1_inf;
er2 = er2_inf;
ur2 = ur2_inf;

% GRID PARAMETERS
ermax = max([er_bg er1 er2 er_sub]);
nmax = sqrt(ermax);
NRESx = 20;
NRESy = 20;
NBUF = [0 0 100 100];
NPML = [0 0 40 40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMAIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOMINAL RESOLUTION
lam0_min = c0/fmax;
dx = lam0_min/nmax/NRESx;
dy = lam0_min/nmax/NRESy;

% SNAP GRID TO CRITICAL DIMENTIONS
Sx = A1+B1+C1;
Nx = 2 * ceil(Sx/dx/2) + 1;
dx = Sx/Nx;
Sy = layer1 + layer2;
Ny = ceil(Sy/dy);
dy = Sy/Ny;

% COMPUTE GRID SIZE
Nx = round(Sx/dx) + 0 + NPML(1) + NPML(2) + NBUF(1) + NBUF(2);
Ny = round(Sy/dy) + 3 + NPML(3) + NPML(4) + NBUF(3) + NBUF(4);

% COMPUTE GRID AXIS
xa = (0:Nx-1)*dx;
ya = (0:Ny-1)*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS TO FREE SPACE
ERxx = er_bg * ones(1,Nx,Ny);
ERyy = er_bg * ones(1,Nx,Ny);
ERzz = er_bg * ones(1,Nx,Ny);
URxx = ur_bg * ones(1,Nx,Ny);
URyy = ur_bg * ones(1,Nx,Ny);
URzz = ur_bg * ones(1,Nx,Ny);

GammaE = zeros(n_res,Nx,Ny);
GammaM = zeros(n_res,Nx,Ny);
WpE = zeros(1,Nx,Ny);
WpM = zeros(1,Nx,Ny);
W0E = zeros(n_res,Nx,Ny);
W0M = zeros(n_res,Nx,Ny);
fE = zeros(n_res,Nx,Ny);
fM = zeros(n_res,Nx,Ny);

% COMPUTE POSITION INDICES
nx1 = 1;
nx2 = nx1 + round((A1+B1+C1)/dx) - 1;
nx_grating1 = 1;
nx_grating2 = nx_grating1 + round(a1/dx) - 1;
nx_grating3 = round(A1/dx);
nx_grating4 = nx_grating3 + 1 + round(b1/dx) -1;
nx_grating5 = round((A1+B1)/dx);
nx_grating6 = nx_grating5 + 1 + round(c1/dx) -1;

ny1 = 2 + NBUF(3) + NPML(3) + 1;
ny2 = ny1 + round(layer1/dy) - 1;
ny3 = ny2 + round(layer2/dy) -1;

% ADD LAYER1 (GRATING)
ERxx(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = er1;
ERyy(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = er1;
ERzz(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = er1;

WpE(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = Wp_E1;

URxx(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = ur1;
URyy(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = ur1;
URzz(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = ur1;

WpM(1,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = Wp_M1;

for i=1:n_res
    
    GammaE(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = Gamma_E1(i);
    W0E(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = W0_E1(i);
    fE(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = f_E1(i);

    GammaM(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = Gamma_M1(i);
    W0M(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = W0_M1(i);
    fM(i,[nx_grating1:nx_grating2,nx_grating3:nx_grating4,nx_grating5:nx_grating6],ny1:ny2) = f_M1(i);

end

% ADD LAYER2
ERxx(1,nx1:nx2,ny2+1:ny3) = er2;
ERyy(1,nx1:nx2,ny2+1:ny3) = er2;
ERzz(1,nx1:nx2,ny2+1:ny3) = er2;

WpE(1,nx1:nx2,ny2+1:ny3) = Wp_E2;

URxx(1,nx1:nx2,ny2+1:ny3) = ur2;
URyy(1,nx1:nx2,ny2+1:ny3) = ur2;
URzz(1,nx1:nx2,ny2+1:ny3) = ur2;

WpM(1,nx1:nx2,ny2+1:ny3) = Wp_M2;

for i=1:n_res
    
    GammaE(i,nx1:nx2,ny2+1:ny3) = Gamma_E2(i);
    W0E(i,nx1:nx2,ny2+1:ny3) = W0_E2(i);
    fE(i,nx1:nx2,ny2+1:ny3) = f_E2(i);

    GammaM(i,nx1:nx2,ny2+1:ny3) = Gamma_M2(i);
    W0M(i,nx1:nx2,ny2+1:ny3) = W0_M2(i);
    fM(i,nx1:nx2,ny2+1:ny3) = f_M2(i);

end

% ADD SUBSTRATE
ERxx(1,nx1:nx2,ny3+1:Ny) = er_sub;
ERyy(1,nx1:nx2,ny3+1:Ny) = er_sub;
ERzz(1,nx1:nx2,ny3+1:Ny) = er_sub;

URxx(1,nx1:nx2,ny3+1:Ny) = ur_sub;
URyy(1,nx1:nx2,ny3+1:Ny) = ur_sub;
URzz(1,nx1:nx2,ny3+1:Ny) = ur_sub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING THE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME STEP
dmin = min([dx dy]);
dt = dmin/(2 * c0);

% SOURCE POSITION
ny_src = NPML(3) + 2;

% COMPUTE SOURCE PARAMETERS
tau = 0.5/fmax;
t0 = 6*tau;

% COMPUTE NUMBER OF TIME STEPS
tprop = nmax*Ny*dy/c0;
t = 2*t0 + 100*tprop;
STEPS = ceil(t/dt);

% COMPUTE THE SOURCE
ta = (0:STEPS-1)*dt;
A = sqrt(ERzz(1,1,ny_src)/URxx(1,1,ny_src));
shift = dy/(2*c0) + dt/2;

Ez_src = exp(-((ta-t0)/tau).^2);
Hx_src = A * exp(-((ta-t0+shift)/tau).^2);
Ex_src = exp(-((ta-t0)/tau).^2);
Hz_src = -A * exp(-((ta-t0+shift)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE THE PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NUMBER OF POINTS ON 2X GRID
Nx2 = 2*Nx;
Ny2 = 2*Ny;

% COMPUTE SIGMA PML PARAMETERS
sigx = zeros(1,Nx2,Ny2);
sigy = zeros(1,Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(1,nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(1,nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end

for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3)- ny +1;
    sigy(1,:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) +ny;
    sigy(1,:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE FIELDS AND RECORD PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISE FOURIER TRANSFORM
K = exp(-1i*2*pi*dt*FREQ);

ErefTE = zeros(Nx,NFREQ);
EtrnTE = zeros(Nx,NFREQ);
SRCTE = zeros(1,NFREQ);

ErefTM = zeros(Nx,NFREQ);
EtrnTM = zeros(Nx,NFREQ);
SRCTM = zeros(1,NFREQ);

% POSITION OF RECORD PLANES
ny_ref = NPML(3) + 1;
ny_trn = Ny - NPML(4);

% COMPUTE REFRACTIVE INDICES IN RECORD PLANES
nref = sqrt(ERzz(1,1,ny_ref)*URxx(1,1,ny_ref));
ntrn = sqrt(ERzz(1,1,ny_trn)*URxx(1,1,ny_trn));

% INITIALISE FIELDS
Ex = zeros(1,Nx,Ny);
Ex_sqz = zeros(Nx,Ny);
Ey = zeros(1,Nx,Ny);
Ez = zeros(1,Nx,Ny);
Ez_sqz = zeros(Nx,Ny);
Dx = zeros(1,Nx,Ny);
Dy = zeros(1,Nx,Ny);
Dz = zeros(1,Nx,Ny);
Hx = zeros(1,Nx,Ny);
Hy = zeros(1,Nx,Ny);
Hz = zeros(1,Nx,Ny);
Bx = zeros(1,Nx,Ny);
By = zeros(1,Nx,Ny);
Bz = zeros(1,Nx,Ny);
JEx = zeros(n_res,Nx,Ny);
JEy = zeros(n_res,Nx,Ny);
JEz = zeros(n_res,Nx,Ny);
JMx = zeros(n_res,Nx,Ny);
JMy = zeros(n_res,Nx,Ny);
JMz = zeros(n_res,Nx,Ny);
Px = zeros(n_res,Nx,Ny);
Py = zeros(n_res,Nx,Ny);
Pz = zeros(n_res,Nx,Ny);
Px_full = zeros(1,Nx,Ny);
Py_full = zeros(1,Nx,Ny);
Pz_full = zeros(1,Nx,Ny);
Mx = zeros(n_res,Nx,Ny);
My = zeros(n_res,Nx,Ny);
Mz = zeros(n_res,Nx,Ny);
Mx_full = zeros(1,Nx,Ny);
My_full = zeros(1,Nx,Ny);
Mz_full = zeros(1,Nx,Ny);

% INITIALIZE CURL
CEx = zeros(1,Nx,Ny);
CEy = zeros(1,Nx,Ny);
CEz = zeros(1,Nx,Ny);
CHx = zeros(1,Nx,Ny);
CHy = zeros(1,Nx,Ny);
CHz = zeros(1,Nx,Ny);

% INITIALIZE INTEGRATION TERMS
ICEx = zeros(1,Nx,Ny);
ICEy = zeros(1,Nx,Ny);
ICEz = zeros(1,Nx,Ny);
IDz = zeros(1,Nx,Ny);
IBz = zeros(1,Nx,Ny);
ICHx = zeros(1,Nx,Ny);
ICHy = zeros(1,Nx,Ny);
ICHz = zeros(1,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE Dx UPDATE COEFFICIENTS
sigDx = sigx(1,2:2:Nx2,1:2:Ny2);
sigDy = sigy(1,2:2:Nx2,1:2:Ny2);
mDx0 = (1/dt) + sigDy/(2*e0);
mDx1 = ((1/dt) - sigDy/(2*e0))./mDx0;
mDx2 = (c0)./mDx0;
mDx3 = (c0*dt*sigDx/e0)./mDx0;

% COMPUTE Dy UPDATE COEFFICIENTS
sigDx = sigx(1,1:2:Nx2,2:2:Ny2);
sigDy = sigy(1,1:2:Nx2,2:2:Ny2);
mDy0 = (1/dt) + sigDx/(2*e0);
mDy1 = ((1/dt) - sigDx/(2*e0))./mDy0;
mDy2 = (c0)./mDy0;
mDy3 = (c0*dt*sigDy/e0)./mDy0;

% COMPUTE Dz UPDATE COEFFICIENTS
sigDx = sigx(1,1:2:Nx2,1:2:Ny2);
sigDy = sigy(1,1:2:Nx2,1:2:Ny2);
sigDz = 0;
mDz0 = (1/dt) + (sigDx+sigDy)/(2*e0) +...
    (sigDx.*sigDy)*dt./4/e0.^2;
mDz1 = ((1/dt) - (sigDx+sigDy)/(2*e0) -...
    (sigDx.*sigDy)*dt./4/e0.^2)./mDz0;
mDz2 = (c0)./mDz0;
mDz3 = (c0*dt*sigDz/e0)./mDz0;
mDz4 = (-sigDx.*sigDy*dt/e0.^2)./mDz0;

% COMPUTE Bx UPDATE COEFFICIENTS
sigBx = sigx(1,1:2:Nx2,2:2:Ny2);
sigBy = sigy(1,1:2:Nx2,2:2:Ny2);
mBx0 = (1/dt) + sigBy/(2*e0);
mBx1 = ((1/dt) - sigBy/(2*e0))./mBx0;
mBx2 = (-N0)./mBx0;
mBx3 = (-N0*dt*sigBx/e0)./mBx0;

% COMPUTE By UPDATE COEFFICIENTS
sigBx = sigx(1,2:2:Nx2,1:2:Ny2);
sigBy = sigy(1,2:2:Nx2,1:2:Ny2);
mBy0 = (1/dt) + sigBx/(2*e0);
mBy1 = ((1/dt) - sigBx/(2*e0))./mBy0;
mBy2 = (-N0)./mBy0;
mBy3 = (-N0*dt*sigBy/e0)./mBy0;

% COMPUTE Bz UPDATE COEFFICIENTS
sigBx = sigx(1,2:2:Nx2,2:2:Ny2);
sigBy = sigy(1,2:2:Nx2,2:2:Ny2);
sigBz = 0;
mBz0 = (1/dt) + (sigBx+sigBy)/(2*e0) +...
    (sigBx.*sigBy)*dt./4/e0.^2;
mBz1 = ((1/dt) - (sigBx+sigBy)/(2*e0) -...
    (sigBx.*sigBy)*dt./4/e0.^2)./mBz0;
mBz2 = (-N0)./mBz0;
mBz3 = (-N0*dt*sigBz/e0)./mBz0;
mBz4 = (-sigBx.*sigBy*dt/e0.^2)./mBz0;

% COMPUTE JEx UPDATE COEFFICIENTS
mJEx0 = (1/dt) + (GammaE/2);
mJEx1 = ((1/dt) - (GammaE/2))./mJEx0;
mJEx2 = (WpE.^2.*fE)./mJEx0;
mJEx3 = (-W0E.^2)./mJEx0;

% COMPUTE JEy UPDATE COEFFICIENTS
mJEy0 = (1/dt) + (GammaE/2);
mJEy1 = ((1/dt) - (GammaE/2))./mJEy0;
mJEy2 = (WpE.^2.*fE)./mJEy0;
mJEy3 = (-W0E.^2)./mJEy0;

% COMPUTE JEz UPDATE COEFFICIENTS
mJEz0 = (1/dt) + (GammaE/2);
mJEz1 = ((1/dt) - (GammaE/2))./mJEz0;
mJEz2 = (WpE.^2.*fE)./mJEz0;
mJEz3 = (-W0E.^2)./mJEz0;

% COMPUTE JMx UPDATE COEFFICIENTS
mJMx0 = (1/dt) + (GammaM/2);
mJMx1 = ((1/dt) - (GammaM/2))./mJMx0;
mJMx2 = (WpM.^2.*fM)./mJMx0;
mJMx3 = (-W0M.^2)./mJMx0;

% COMPUTE JMy UPDATE COEFFICIENTS
mJMy0 = (1/dt) + (GammaM/2);
mJMy1 = ((1/dt) - (GammaM/2))./mJMy0;
mJMy2 = (WpM.^2.*fM)./mJMy0;
mJMy3 = (-W0M.^2)./mJMy0;

% COMPUTE JMz UPDATE COEFFICIENTS
mJMz0 = (1/dt) + (GammaM/2);
mJMz1 = ((1/dt) - (GammaM/2))./mJMz0;
mJMz2 = (WpM.^2.*fM)./mJMz0;
mJMz3 = (-W0M.^2)./mJMz0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDTD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAIN LOOP - ITERATIVE OVER TIME

for T = 1 : STEPS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TE mode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Mx_full = zeros(1,Nx,Ny);    
    My_full = zeros(1,Nx,Ny);
    
    Pz_full = zeros(1,Nx,Ny);
            
    for i = 1 : n_res
    
        % UPDATE JMx,JMy FROM Mx,My AND Hx,Hy 
        JMx(i,:,:) = mJMx1(i,:,:).*JMx(i,:,:) +...
            mJMx2(i,:,:).*Hx + mJMx3(i,:,:).*Mx(i,:,:);
        JMy(i,:,:) = mJMy1(i,:,:).*JMy(i,:,:) +...
            mJMy2(i,:,:).*Hy + mJMy3(i,:,:).*My(i,:,:);
            
        % UPDATE Mx,My FROM JMx,JMy
        Mx(i,:,:) = Mx(i,:,:) + dt.*JMx(i,:,:);
        My(i,:,:) = My(i,:,:) + dt.*JMy(i,:,:);
                
        Mx_full = Mx_full + Mx(i,:,:);
        My_full = My_full + My(i,:,:);
            
    end 
    
    % COMPUTE CEx
    for ny = 1 :Ny-1
        for nx = 1 : Nx
            CEx(1,nx,ny) = (Ez(1,nx,ny+1)-Ez(1,nx,ny))/dy;
        end
    end
    for nx = 1:Nx
        CEx(1,nx,Ny) = (Ez(1,nx,1)-Ez(1,nx,Ny))/dy;
    end
   
    % COMPUTE CEy
    for nx = 1 :Nx-1
        for ny = 1 :Ny
            CEy(1,nx,ny) =-(Ez(1,nx+1,ny)-Ez(1,nx,ny))/dx;
        end
    end
    for ny = 1:Ny
        CEy(1,Nx,ny) = -(Ez(1,1,ny)-Ez(1,Nx,ny))/dx;
    end
        
    % TF/SF
    CEx(1,:,ny_src-1) = CEx(1,:,ny_src-1) - Ez_src(T)/dy;
     
    % UPDATE ICEx,ICEy
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    
    % UPDTAE B FROM E
    Bx = mBx1.*Bx + mBx2.*CEx + mBx3.*ICEx;
    By = mBy1.*By + mBy2.*CEy + mBy3.*ICEy; 
    
    % UPDATE H FROM B AND M
    Hx = (Bx-Mx_full)./URxx/u0;
    Hy = (By-My_full)./URyy/u0;   
    
    % COMPUTE CHz
    CHz(1,1,1) = (Hy(1,1,1) - Hy(1,Nx,1))/dx - (Hx(1,1,1) -Hx(1,1,Ny))/dy;
    for nx = 2 : Nx
        CHz(1,nx,1) = (Hy(1,nx,1)-Hy(1,nx-1,1))/dx - (Hx(1,nx,1)-Hx(1,nx,Ny))/dy;
    end
    for ny = 2 : Ny
        CHz(1,1,ny) = (Hy(1,1,ny) - Hy(1,Nx,ny))/dx - ...
            (Hx(1,1,ny) - Hx(1,1,ny-1))/dy;
        for nx = 2 : Nx
            CHz(1,nx,ny) = (Hy(1,nx,ny) - Hy(1,nx-1,ny))/dx - ...
                (Hx(1,nx,ny) - Hx(1,nx,ny-1))/dy;
        end
    end
    
    % TF/SF
    CHz(1,:,ny_src) = CHz(1,:,ny_src) + Hx_src(T)/dy;
    
    for i = 1 : n_res
           
        % UPDATE JEz FROM Pz AND Ez
        JEz(i,:,:) = mJEz1(i,:,:).*JEz(i,:,:) +...
            mJEz2(i,:,:).*Ez + mJEz3(i,:,:).*Pz(i,:,:);
    
        % UPDATE Pz FROM JEz
        Pz(i,:,:) = Pz(i,:,:) + dt.*JEz(i,:,:);
            
        Pz_full = Pz_full + Pz(i,:,:);
            
    end
        
    % UPDATE IDz
    IDz = IDz + Dz;
    
    % UPDATE ICHz
    ICHz = ICHz + CHz;
    
    % UPDTAE Dz FROM Hx,Hy
    Dz = mDz1.*Dz + mDz2.*CHz + mDz3.*ICHz + mDz4.*IDz;
    
    % UPDATE E FROM D AND P
    Ez = (Dz - Pz_full)./ERzz;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TM MODE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Mz_full = zeros(1,Nx,Ny);    
        
    Px_full = zeros(1,Nx,Ny);
    Py_full = zeros(1,Nx,Ny);
    
    for i = 1 : n_res
    
        % UPDATE JMz FROM Mz AND Hz 
        JMz(i,:,:) = mJMz1(i,:,:).*JMz(i,:,:) +...
            mJMz2(i,:,:).*Hz + mJMz3(i,:,:).*Mz(i,:,:);
                    
        % UPDATE Mz FROM JMz
        Mz(i,:,:) = Mz(i,:,:) + dt.*JMz(i,:,:);
                
        Mz_full = Mz_full + Mz(i,:,:);
                    
    end 
    
    % COMPUTE CEz
    for ny = 1 : Ny-1
        for nx = 1 : Nx-1
            CEz(1,nx,ny) = (Ey(1,nx+1,ny) - Ey(1,nx,ny))/dx - ...
                (Ex(1,nx,ny+1) - Ex(1,nx,ny))/dy;
        end
        CEz(1,Nx,ny) = (Ey(1,1,ny) - Ey(1,Nx,ny))/dx - ...
            (Ex(1,Nx,ny+1) - Ex(1,Nx,ny))/dy;
    end
    for nx = 1 : Nx-1
        CEz(1,nx,Ny) = (Ey(1,nx+1,Ny)-Ey(1,nx,Ny))/dx - (Ex(1,nx,1)-Ex(1,nx,Ny))/dy;
    end
    CEz(1,Nx,Ny) = (Ey(1,1,Ny) - Ey(1,Nx,Ny))/dx - (Ex(1,Nx,1) -Ex(1,Nx,Ny))/dy;
    
    % TF/SF
    CEz(1,:,ny_src-1) = CEz(1,:,ny_src-1) + Ex_src(T)/dy;
        
    % UPDATE ICEz
    ICEz = ICEz + CEz;

    % UPDATE IBz
    IBz = IBz + Bz;
    
    % UPDTAE Bz FROM Ex,Ey
    Bz = mBz1.*Bz + mBz2.*CEz + mBz3.*ICEz + mBz4.*IBz;
     
    % UPDATE Hz FROM Bz AND Mz
    Hz = (Bz-Mz_full)./URzz/u0;
       
    % COMPUTE CHx
    for nx = 1:Nx
        CHx(1,nx,1) = (Hz(1,nx,1)-Hz(1,nx,Ny))/dy;
    end
    for ny = 2 :Ny
        for nx = 1 : Nx
            CHx(1,nx,ny) = (Hz(1,nx,ny)-Hz(1,nx,ny-1))/dy;
        end
    end
    
    % COMPUTE CHy
    for ny = 1:Ny
        CHy(1,1,ny) = -(Hz(1,1,ny)-Hz(1,Nx,ny))/dx;
    end 
    for nx = 2 :Nx
        for ny = 1 :Ny
            CHy(1,nx,ny) =-(Hz(1,nx,ny)-Hz(1,nx-1,ny))/dx;
        end
    end
       
    % TF/SF
    CHx(1,:,ny_src) = CHx(1,:,ny_src) - Hz_src(T)/dy;
    
    for i = 1 : n_res
           
        % UPDATE JEx,JEy FROM Px,Py AND Ex,Ey
        JEx(i,:,:) = mJEx1(i,:,:).*JEx(i,:,:) +...
            mJEx2(i,:,:).*Ex + mJEx3(i,:,:).*Px(i,:,:);
        JEy(i,:,:) = mJEy1(i,:,:).*JEy(i,:,:) +...
            mJEy2(i,:,:).*Ey + mJEy3(i,:,:).*Py(i,:,:);
       
        % UPDATE Px,Py FROM JEx,JEy
        Px(i,:,:) = Px(i,:,:) + dt.*JEx(i,:,:);
        Py(i,:,:) = Py(i,:,:) + dt.*JEy(i,:,:);
    
        Px_full = Px_full + Px(i,:,:);
        Py_full = Py_full + Py(i,:,:);
            
    end
        
    % UPDATE ICHx,ICHy
    ICHx = ICHx + CHx;
    ICHy = ICHy + CHy;
    
    % UPDTAE Dx,Dy FROM Hz
    Dx = mDx1.*Dx + mDx2.*CHx + mDx3.*ICHx;
    Dy = mDy1.*Dy + mDy2.*CHy + mDy3.*ICHy;

    % UPDATE Ex,Ey FROM Dx,Dy AND Px,Py
    Ex = (Dx - Px_full)./ERxx;
    Ey = (Dy - Py_full)./ERyy;
           
     % UPDATE FOURIER TRANSFORMS
    Ex_sqz = squeeze(Ex);
    Ez_sqz = squeeze(Ez);
    
    for nfreq = 1 : NFREQ 
       ErefTE(:,nfreq) = ErefTE(:,nfreq) + (K(nfreq)^T)*Ez_sqz(:,ny_ref)*dt;
       EtrnTE(:,nfreq) = EtrnTE(:,nfreq) + (K(nfreq)^T)*Ez_sqz(:,ny_trn)*dt;
       SRCTE(1,nfreq) = SRCTE(1,nfreq) + (K(nfreq)^T)*Ez_src(T)*dt;
    
       ErefTM(:,nfreq) = ErefTM(:,nfreq) + (K(nfreq)^T)*Ex_sqz(:,ny_ref)*dt;
       EtrnTM(:,nfreq) = EtrnTM(:,nfreq) + (K(nfreq)^T)*Ex_sqz(:,ny_trn)*dt;
       SRCTM(1,nfreq) = SRCTM(1,nfreq) + (K(nfreq)^T)*Ex_src(T)*dt;
    end
 
    % SHOW STATUS
    if ~mod(T, update_plot)
        
        % SHOW FIELDS
        subplot(2,5,1);
        draw2D_BS(xa,ya,squeeze(ERzz),squeeze(Ez),NPML,CMAP);
        axis equal tight;
        title([num2str(T) ' of ' num2str(STEPS)]);
        ylabel('TE mode');
        
        subplot(2,5,6);
        draw2D_BS(xa,ya,squeeze(ERxx),squeeze(Ex),NPML,CMAP);
        axis equal tight;
        ylabel('TM mode');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% REFLECTANCE AND TRANSMITTANCE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
        % INITIALIZE REFLECTANCE AND TRANSMITTANCE
        REFTE = zeros(1,NFREQ);
        TRNTE = zeros(1,NFREQ);

        REFTM = zeros(1,NFREQ);
        TRNTM = zeros(1,NFREQ);

        % LOOP OVER FREQUENCY
        for nfreq = 1 : NFREQ
            
            % COMPUTE WAVE VECTOR COMPONENTS
            lam0 = c0/FREQ(nfreq);
            k0 = 2*pi/lam0;
            kxinc = 0;
            kyinc = k0*nref;
            m = (-floor(Nx/2):floor(Nx/2))';
            kx = -2*pi*m/Sx;
            kyR = sqrt((k0*nref)^2 - kx.^2);
            kyT = sqrt((k0*ntrn)^2 - kx.^2);

            % COMPUTE REFLECTANCE
            refTE = ErefTE(:,nfreq)/SRCTE(1,nfreq);
            refTE = fftshift(fft(refTE))/Nx;
            refTE = real(kyR/kyinc).*abs(refTE).^2;
            REFTE(nfreq) = sum(refTE);

            refTM = ErefTM(:,nfreq)/SRCTM(1,nfreq);
            refTM = fftshift(fft(refTM))/Nx;
            refTM = real(kyR/kyinc).*abs(refTM).^2;
            REFTM(nfreq) = sum(refTM);

            % COMPUTE TRANSMITTANCE
            trnTE = EtrnTE(:,nfreq)/SRCTE(1,nfreq);
            trnTE = fftshift(fft(trnTE))/Nx;
            trnTE = real(kyT/kyinc).*abs(trnTE).^2;
            TRNTE(nfreq) = sum(trnTE);

            trnTM = EtrnTM(:,nfreq)/SRCTM(1,nfreq);
            trnTM = fftshift(fft(trnTM))/Nx;
            trnTM = real(kyT/kyinc).*abs(trnTM).^2;
            TRNTM(nfreq) = sum(trnTM);

        end
        
        % COMPUTE ENERGY CONSERVATION
        CONTE = REFTE + TRNTE;

        CONTM = REFTM + TRNTM;
            
        % PLOT REFLECTION/TRANSMISSION SPECTRUM
        subplot(2,5,(2:5));         
        plot(FREQ,REFTE,'-R','LineWidth',2);
        hold on;
        plot(FREQ,TRNTE,'-B','LineWidth',2);
        plot(FREQ,CONTE,'-K','LineWidth',1);
        axis([FREQ(1) FREQ(NFREQ) -0.1 1.1]);
        xlabel('Frequency (Hz)'); ylabel('Normalized intensity (a.u)');
        title('TE mode');
        legend('Reflectance','Transmittance','Conservation');
        hold off;
                    
        subplot(2,5,(7:10)); 
        plot(FREQ,REFTM,'-R','LineWidth',2);
        hold on;
        plot(FREQ,TRNTM,'-B','LineWidth',2);
        plot(FREQ,CONTM,'-K','LineWidth',1);
        axis([FREQ(1) FREQ(NFREQ) -0.1 1.1]);
        xlabel('Frequency (Hz)'); ylabel('Normalized intensity (a.u)');
        title('TM mode');
        legend('Reflectance','Transmittance','Conservation');
        hold off;
        
        drawnow;
       
    end
    
    drawnow;
    
end
        
%% GENERATE A PROFESSIONAL LOOKING PLOT

fig2 = figure('numbertitle','off','name','Spectrum',...
    'color','white','position',scrnsize);
h = plot(c0./FREQ./nm,100*REFTE,'-r','LineWidth',2);
hold on;
plot(c0./FREQ./nm,100*REFTM,'--r','LineWidth',2);
plot(c0./FREQ./nm,100*TRNTE,'-b','LineWidth',2);
plot(c0./FREQ./nm,100*TRNTM,'--b','LineWidth',2);
plot(c0./FREQ./nm,100*CONTE,'-k','LineWidth',1);
plot(c0./FREQ./nm,100*CONTM,'--k','LineWidth',1);
hold off;
axis([c0./FREQ(NFREQ)./nm c0./FREQ(1)./nm 0 105]);
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
h = legend('Reflectance - TE','Reflectance - TM',...
    'Transmittance - TE','Transmittance - TM','Conservation - TE','Conservation - TM');
set(h,'Location','NorthEast');
xlabel('Wavelength/nm');
ylabel('%','Rotation',0,'HorizontalAlignment','right');


fig3 = figure('numbertitle','off','name','Spectrum',...
    'color','white','position',scrnsize);
h = plot(c0./FREQ./nm,100*TRNTE,'-r','LineWidth',2);
hold on;
plot(c0./FREQ./nm,100*TRNTM,'--b','LineWidth',2);
hold off;
axis([c0./FREQ(NFREQ)./nm c0./FREQ(1)./nm 0 105]);
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
h = legend(...
    'Transmittance - TE','Transmittance - TM');
set(h,'Location','NorthEast');
xlabel('Wavelength/nm');
ylabel('%','Rotation',0,'HorizontalAlignment','right');