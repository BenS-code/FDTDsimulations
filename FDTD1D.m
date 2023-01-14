%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINITE DIFFERENCE TIME DOMAIN SIMULATION - 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
The 1D simulation (upgraded from Raymond Rampf's lectures) solves the
reflection & transmission spectrum of a pulse
impinges perpendicularly on a slab with material properties as summarized
in the pdf file in folder "FDTD_equations".

Spatial grid is based on "YEE grid" 

Material properties, frequency range and the number of layers can be
modified.
%}


%% INITIALIZE MATLAB

close all; clear all; clc;

%% INITIALIZE FIGURE 

scrnsize = get(0,'screensize');
fig1 = figure('numbertitle','off','name','FDTD Simulation',...
    'color','white','position',scrnsize);

%% UNITS

m = 1;          % distance units
cm = 1e-2 * m;
mm = 1e-3 * m;
um = 1e-6 * m;
nm = 1e-9 * m;

s = 1;          % time units
ms = 1e-3 * s;
us = 1e-6 * s;
ns = 1e-9 * s;
ps = 1e-12 * s;

Hz = 1;         % frequency units
kHz = 1e3 * Hz;
MHz = 1e6 * Hz;
GHz = 1e9 * Hz;
THz = 1e12 * Hz;
PHz = 1e15 * Hz;
             
H = 1;           % Henry units
Ohm = 1;         % Resistance units

%% CONSTANTS

e0 = 8.8541878176 * 1e-12 * Hz/m;   % permitivity
u0 = 1.2566370614 * 1e-6 * H/m;     % permiability
c0 = 1/sqrt(e0*u0) * m/s;           % speed of light in vacuum
N0 = sqrt(u0/e0) * Ohm;             % impedance of vacuum

%% SIMULATION PROPERTIES

% FREQUENCY RANGE
fmin = 0.1 * PHz;                   % minimum frequency
fmax = 1.0 * PHz;                   % maximum frequency
NFREQ = 500;                        % number of frequencies
FREQ = linspace(fmin,fmax,NFREQ);   % frequency array
        
% MATERIAL PROPERTIES 
% dispersive electric material (see FDTD_equantions):
Wp_E = 13.72e15;                            
Gamma_E = [0.08052 0.3661 0.5241].*1e15;    
W0_E = [0 0.6305 1.261].*1e15;
f_E = [0.76 0.024 0.01];

er_inf = 8.9;

% dispersive magnetic material (see FDTD_equantions): 
Wp_M = 0;
Gamma_M = [0 0 0];
W0_M = [0 0 0];
f_M = [0 0 0];
u_inf = 0;

n_res = length(f_E);                % number of resonances in the material

% SLAB PROPERTIES
dslab = 100 * nm;
erair= 1.0;
ersubstrate = 1.0;
erslab = er_inf;
urslab = u_inf;

% GRID PARAMETERS
ermax = max([erair erslab ersubstrate]);
nmax = sqrt(ermax);
NLAM = 40;
NDIM = 4;
NBUFZ = [200 200];
NPML = [20 20];

update_plot = 50;        % update plot every "update_plot" time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMAIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOMINAL RESOLUTION
lam0_min = c0/fmax;
dz1 = lam0_min/nmax/NLAM;
dz2 = dslab/NDIM;
dz = min([dz1 dz2]); 

% SNAP GRID TO CRITICAL DIMENTIONS
nz = ceil(dslab/dz);
dz = dslab/nz;

% COMPUTE GRID SIZE
Nz = round(dslab/dz) + sum(NBUFZ) + 3 + sum(NPML);

% COMPUTE GRID AXIS
za = (0:Nz-1)*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS TO FREE SPACE
ER = erair * ones(1,Nz);
UR = ones(1,Nz);
GammaE = zeros(n_res,Nz);
GammaM = zeros(n_res,Nz);
WpE = zeros(1,Nz);
WpM = zeros(1,Nz);
W0E = zeros(n_res,Nz);
W0M = zeros(n_res,Nz);
fE = zeros(n_res,Nz);
fM = zeros(n_res,Nz);

% COMPUTE POSITION INDICES
nz1 = 2 + NBUFZ(1) + NPML(1) + 1;
nz2 = nz1 + round(dslab/dz) - 1;

% ADD THE SLAB 
ER(nz1:nz2) = erslab;
WpE(nz1:nz2) = Wp_E;

UR(nz1:nz2) = urslab;
WpM(nz1:nz2) = Wp_M;

for i=1:n_res
    
    GammaE(i,nz1:nz2) = Gamma_E(i);
    W0E(i,nz1:nz2) = W0_E(i);
    fE(i,nz1:nz2) = f_E(i);

    GammaM(i,nz1:nz2) = Gamma_M(i);
    W0M(i,nz1:nz2) = W0_M(i);
    fM(i,nz1:nz2) = f_M(i);

end

% ADD SUBSTRATE
ER(nz2+1:Nz) = ersubstrate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTING THE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute time steps (dt)

nbc = sqrt(UR(NPML(1)+1)*ER(NPML(1)+1));
dt = nbc*dz/(2*c0);

% COMPUTE SOURCE PARAMETERS
tau = 0.5/fmax;
t0 = 6*tau;

% COMPUTE NUMBER OF TIME STEPS
tprop = nmax*Nz*dz/c0;
t = 2*t0 + 5*tprop;
STEPS = ceil(t/dt);

% COMPUTE THE SOURCE
t = (0:STEPS-1)*dt;
shift = dz/(2*c0) +dt/2;
nz_src = 2 + NPML(1);
Ey_src = exp(-((t-t0)/tau).^2);
A = -sqrt(ER(nz_src)/UR(nz_src));
Hx_src = A*exp(-((t-t0+shift)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE THE PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NUMBER OF POINTS ON 2X GRID
Nz2 = 2*Nz;

% COMPUTE SIGMA PML PARAMETERS
sigz = zeros(1,Nz2);
for nz = 1 : 2*NPML(1)
    nz1 = 2*NPML(1)- nz +1;
    sigz(1,nz1) = (0.5*e0/dt)*(nz/2/NPML(1))^3;
end
for nz = 1 : 2*NPML(2)
    nz1 = Nz2 - 2*NPML(2) + nz;
    sigz(1,nz1) = (0.5*e0/dt)*(nz/2/NPML(2))^3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE FDTD PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POSITION OF RECORD PLANES
n_ref = NPML(1)+1;
n_trn = Nz-NPML(2);

% INITIALISE FOURIER TRANSFORM
K = exp(-1i*2*pi*dt*FREQ);
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMUTE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE Bx UPDTAE COEFFICIENTS
sigBz = sigz(2:2:Nz2);
mBx0 = (1/dt) + sigBz/(2*e0);
mBx1 = ((1/dt) - sigBz/(2*e0))./mBx0;
mBx2 = (-N0)./mBx0;

sigDz = sigz(1:2:Nz2);
mDy0 = (1/dt) + sigDz/(2*e0);
mDy1 = ((1/dt) - sigDz/(2*e0))./mDy0;
mDy2 = (c0)./mDy0;

mJEy0 = (1/dt) + (GammaE/2);
mJEy1 = ((1/dt) - (GammaE/2))./mJEy0;
mJEy2 = (WpE.^2.*fE)./mJEy0;
mJEy3 = (-W0E.^2)./mJEy0;

mJMx0 = (1/dt) + (GammaM/2);
mJMx1 = ((1/dt) - (GammaM/2))./mJMx0;
mJMx2 = (WpM.^2.*fM)./mJMx0;
mJMx3 = (-W0M.^2)./mJMx0;

% INITIALISE FIELDS
Ey = zeros(1,Nz);
Dy = zeros(1,Nz);
JEy = zeros(n_res,Nz);
Py = zeros(n_res,Nz);
Py_full = zeros(1,Nz);
Hx = zeros(1,Nz);
Bx = zeros(1,Nz);
JMx = zeros(n_res,Nz);
Mx = zeros(n_res,Nz);
Mx_full = zeros(1,Nz);

CEx = zeros(1,Nz);
CHy = zeros(1,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDTD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MAIN LOOP - ITERATIVE OVER TIME

for T = 1 : STEPS
     
    Mx_full = zeros(1,Nz);
    Py_full = zeros(1,Nz);
    for i = 1 : n_res
    
        % UPDATE JM FROM M AND H 
        JMx(i,:) = mJMx1(i,:).*JMx(i,:) +...
            mJMx2(i,:).*Hx + mJMx3(i,:).*Mx(i,:);
    
        % UPDATE M FROM JM
        Mx(i,:) = Mx(i,:) + dt.*JMx(i,:);
        
        Mx_full = Mx_full + Mx(i,:);
    
    end 
    
        % COMPUTE CEx
        for nz = 1 : Nz-1
            CEx(nz) = -(Ey(nz+1)-Ey(nz))/dz;
        end
        CEx(Nz) = -(0-Ey(Nz))/dz;
    
        % TF/SF
        CEx(nz_src-1) = CEx(nz_src-1) + Ey_src(T)/dz;
    
        % UPDTAE B FROM E
        Bx = mBx1.*Bx + mBx2.*CEx;
        
        % UPDATE H FROM B AND M
        Hx = (Bx-Mx_full)./UR/u0;
        
        % COMPUTE CHy
        CHy(1) = (Hx(1)-0)/dz;
        for nz = 2 : Nz
            CHy(nz) = (Hx(nz)-Hx(nz-1))/dz;
        end
    
        % TF/SF
        CHy(nz_src) = CHy(nz_src) - Hx_src(T)/dz;
        
    for i = 1 : n_res
           
        % UPDATE JE FROM P AND E
        JEy(i,:) = mJEy1(i,:).*JEy(i,:) +...
            mJEy2(i,:).*Ey + mJEy3(i,:).*Py(i,:);
    
        % UPDATE P FROM JE
        Py(i,:) = Py(i,:) + dt.*JEy(i,:);
            
        Py_full = Py_full + Py(i,:);
            
    end
        
    % UPDTAE D FROM H
    Dy = mDy1.*Dy + mDy2.*CHy;
    
    % UPDATE E FROM D AND P
    Ey = (Dy - Py_full)./ER;
    
    % UPDATE FOURIER TRANSFORMS
    for nf = 1 : NFREQ
    REF(nf) = REF(nf) + (K(nf)^T)*Ey(n_ref);
    TRN(nf) = TRN(nf) + (K(nf)^T)*Ey(n_trn);
    SRC(nf) = SRC(nf) + (K(nf)^T)*Ey_src(T);
    end
    
    % SHOW STATUS
    if ~mod(T,update_plot)
        
        % SHOW FIELDS
        subplot(211);
        draw1D_BS(ER,Ey,Hx,dz,NPML);
        xlim([dz Nz*dz]);
        xlabel('z');
        title(['Fields at step ' num2str(T) ' of ' num2str(STEPS)]);
        
        % TEMPORAY NORMALIZATION
        R = abs(REF./SRC).^2;
        Tr = abs(TRN./SRC).^2;
        
        % SHOW REF AND TRN
        subplot(212);
        plot(FREQ,R,'-r');
        hold on;
        plot(FREQ,Tr,'-B');
        plot(FREQ,R+Tr,':K');
        xlim([FREQ(1) FREQ(NFREQ)]);
        ylim([-0.1 1.1]);
        xlabel('Frequency (Hz)');
        title('Reflactance and Transmittance');
        hold off;
        % DRAW GRAPHICS
        drawnow;
    end
    drawnow;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE A PROFESSIONAL LOOKING PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig2 = figure('numbertitle','off','name','Spectrum',...
    'color','white','position',scrnsize);
h = plot (c0./FREQ./nm,100*R,'-r','LineWidth',2);
hold on;
plot(c0./FREQ./nm,100*Tr,'-b','LineWidth',2);
plot(c0./FREQ./nm,100*(R+Tr),':k','LineWidth',2);
hold off;
axis([c0./FREQ(NFREQ)./nm c0./FREQ(1)./nm 0 105]);
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
h = legend('Reflectance','Transmittance','Conservation');
set(h,'Location','NorthEast');
xlabel('Wavelength/nm');
ylabel('%','Rotation',0,'HorizontalAlignment','right');