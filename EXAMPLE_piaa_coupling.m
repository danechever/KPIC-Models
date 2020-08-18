%% EXAMPLE_piaa_coupling.m
% Script to simulate the image in the KPIC mask testbed. 
% Uses functions from VFNlib https://github.com/danechever/VFN-Simulations

clear; 

addpath('util')

% Add the VFNlib and falco to your path 
path2vfnlib = '~/Desktop/VFN/VFN-Simulations/VFNlib';
addpath(path2vfnlib);


%% Configuration 

%%- Choose mask mode ('piaa', 'apod', 'vfn', or 'none')
mode = 'piaa';
% mode = 'apod';
% mode = 'vfn';
% mode = 'none';

%%- Fiber properties (properties store in fiber_props struct)
% Parameters for Thorlabs SM2000
% link: https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=949&pn=SM2000#4668
source_fiber_props.core_rad = 11e-6/2;% Core radius [meters]
source_fiber_props.n_core = 1.4436;% core index
source_fiber_props.n_clad = 1.4381;% cladding index

% Parameters for the ZBLAN fiber in the FIU
fiu_fiber_props.core_rad = 6.5e-6;% Core radius [meters]
fiu_fiber_props.NA = 0.175;% Measured to be 0.175+-.01
fiu_fiber_props.type = 'bessel'; % This one is more accurate 

%%- Optical system (properties store in sys_props struct)
% The term "mask" refers to the piaa, apod, or vfn, depending on 'mode'
% For the piaa, the mask position is the position of the first aspheric
% surface (i.e. the backside of the first lens)
sys_props.z_FL = 200e-3;% distance from fiber to first lens [meters]
sys_props.z_LI = 50e-3;% distance from the lens to the iris [meters]
sys_props.z_IM = 10e-3;% distance from the iris to the "mask" (or first piaa surface) [meters]
% Note: it is assumed that the distance between the iris and the second lens is
% one focal length.
sys_props.f_lens1 = 200e-3; % focal length of the first lens [meters]
sys_props.f_lens2 = 1.0; % focal length of the second lens [meters]
sys_props.f_fiber = 36.6e-3; % focal length of the lens for fiber coupling [meters]
sys_props.D_I = 12.3e-3;% diameter of the iris [meters]
sys_props.pupil_shape = 'kecklab'; %Pupil shape ('circ' or 'keck' or 'kecklab')

% Transfer matrix between fiber and iris (propagate, lens, propagate)
TM = @(z1,z2,f) [1, z2; 0  , 1]*... %second propagation
                [1, 0; -1/f, 1]*... %lens
                [1, z1; 0  , 1];    %first propagation
            
% Store transfer matrix in sys_props struct (wavelength independent)
sys_props.TM = TM(sys_props.z_FL,sys_props.z_LI,sys_props.f_lens1);

%%- PIAA properties 
piaa_props.filename = 'piaa/KPIAA-01.csv';% file containing sag profiles 
piaa_props.L = 50e-3; % Distance between PIAA surfaces [meters]
piaa_props.material = 'CaF2'; % piaa lens material 
piaa_props.Dstop = 20e-3; % Diameter of the output beam [meters]

% %%- APOD properties (can load from fits file or coefficients) 
% apod_props.D = sys_props.D_I;% diameter of the iris [meters]
% apod_props.type = 'fits'; % Load from fits file. % Requires dx = 25um (i.e. N_iris = D_I/25e-6)
% % apod_props.type = 'coeffs'; % Load analytical expression and coefficients. 
% apod_props.whichdesign = 'KAPOD-02'; % KAPOD-01 was faulty. KAPOD-02 is the replacement.
% 
% %%- VFN properties
% vfn_props.charge = 2; 

%%- Other beam properties
N_iris = 1000; % Number of samples across the beam at the iris
% Wavelengths stored as a vector (e.g. wvls = [2.0e-6 2.2e-6 2.4e-6];)
wvls = (2:0.2:2.4)*1e-6;% wavelengths [meters]
%wvls = (3.4:0.2:4.0)*1e-6;% wavelengths [meters]
% Beam sampling at camera
N_lambdaFnum = 4; % Number of samples per lambda*F# at the camera 
N_pix = 64;  % Number of pixels in final image
dx_pix = 45e-6; % Sample spacing in final image, or pixel pitch [meters]

outfile_label = [sys_props.pupil_shape,'_',mode,'_KPIAA-01_1'];

if(strcmpi(mode,'apod') && strcmpi(apod_props.type,'fits'))
    N_iris = apod_props.D/25e-6;
    disp('Overriding N_iris to match apodizer file.');
end
%% Get E-field at the iris 

N_grid = N_iris*N_lambdaFnum;
[E_iris,coords,w_I] = getEinput(source_fiber_props,sys_props,wvls,N_iris,N_grid);

%%- Plot the E-field at the iris 

xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

% figure(1);
% for wvl_index = 1:numel(wvls)
%     
%     % amplitude 
%     ax1=subplot(numel(wvls),2,2*wvl_index-1);
%     imagesc(xvals/N_iris,yvals/N_iris,abs(E_iris(:,:,wvl_index)).^2);
%     axis image;set(gca,'ydir','normal');
%     axis([-0.5 0.5 -0.5 0.5]);
%     colorbar;colormap(ax1,gray);
%     
%     % phase 
%     ax2=subplot(numel(wvls),2,2*wvl_index);
%     imagesc(xvals/N_iris,yvals/N_iris,angle(E_iris(:,:,wvl_index)));
%     axis image;set(gca,'ydir','normal');
%     axis([-0.5 0.5 -0.5 0.5]);
%     caxis([-pi pi]);colorbar;colormap(ax2,hsv);
%     
% end

%% Get E-field at the mask 

dx = sys_props.D_I/N_iris; % sample spacing [meters]
E_b4mask = propagateFresnel(E_iris,sys_props.z_IM,wvls,dx);

%%- Plot the E-field just before the mask

% figure(1);
% for wvl_index = 1:numel(wvls)
%     
%     % amplitude 
%     ax1=subplot(numel(wvls),2,2*wvl_index-1);
%     imagesc(xvals/N_iris,yvals/N_iris,abs(E_b4mask(:,:,wvl_index)).^2);
%     axis image;set(gca,'ydir','normal');
%     axis([-0.5 0.5 -0.5 0.5]);
%     colorbar;colormap(ax1,gray);
%     
%     % phase 
%     ax2=subplot(numel(wvls),2,2*wvl_index);
%     imagesc(xvals/N_iris,yvals/N_iris,angle(E_b4mask(:,:,wvl_index)));
%     axis image;set(gca,'ydir','normal');
%     axis([-0.5 0.5 -0.5 0.5]);
%     caxis([-pi pi]);colorbar;colormap(ax2,hsv);
%     
% end

%% Apply mask 

%%- Apply the mask corresponding to the selected mode. 
switch lower(mode)
    case 'piaa'
        addpath('piaa')
        E_mask = applyPIAA(E_b4mask,piaa_props,wvls,dx);
    case 'apod'
        addpath('apod');
        E_mask = applyAPOD(E_b4mask,apod_props,dx);
    case 'vfn'
        addpath('vfn');
        E_mask = applyVORTEX(E_b4mask,vfn_props);
    case 'none'
        E_mask = E_b4mask;
    otherwise
        error('mode not recognized. options are piaa, apod, vfn, or none.')
end

% Calculate the effective (virtual) E-field at the iris 
E_eff = propagateFresnel(E_mask,-sys_props.z_IM,wvls,dx);

figure(1);
for wvl_index = 1:numel(wvls)
    
    % amplitude 
    ax1=subplot(numel(wvls),2,2*wvl_index-1);
    imagesc(xvals/N_iris,yvals/N_iris,abs(E_eff(:,:,wvl_index)).^2);
    axis image;set(gca,'ydir','normal');
    axis([-0.5 0.5 -0.5 0.5]);
    colorbar;colormap(ax1,gray);
    
    % phase 
    ax2=subplot(numel(wvls),2,2*wvl_index);
    imagesc(xvals/N_iris,yvals/N_iris,angle(E_eff(:,:,wvl_index)));
    axis image;set(gca,'ydir','normal');
    axis([-0.5 0.5 -0.5 0.5]);
    caxis([-pi pi]);colorbar;colormap(ax2,hsv);
    
end

%% FT to image plane, get camera image

[E_fp ,I_fp ] = getEfp(E_eff ,E_iris,sys_props,wvls,dx,dx_pix,N_pix);
[E_fp0,I_fp0] = getEfp(E_iris,E_iris,sys_props,wvls,dx,dx_pix,N_pix);

xvals_cam = -N_pix/2:N_pix/2-1;yvals_cam = xvals_cam';
figure(901);
for wvl_index = 1:numel(wvls)
    
    % amplitude 
    ax1=subplot(numel(wvls),2,2*wvl_index-1);
    imagesc(xvals_cam/N_lambdaFnum,yvals_cam/N_lambdaFnum,log10(I_fp0(:,:,wvl_index)));
    axis image;set(gca,'ydir','normal');
    colorbar;colormap(hot);
    caxis([-4 0]);
    title('PSF without mask')
    
    % phase 
    ax2=subplot(numel(wvls),2,2*wvl_index);
    imagesc(xvals_cam/N_lambdaFnum,yvals_cam/N_lambdaFnum,log10(I_fp(:,:,wvl_index)));
    axis image;set(gca,'ydir','normal');
    colorbar;colormap(hot);
    caxis([-4 0]);
    title('PSF with mask')
    
end

%% Get coupling efficiencies 

etas = getCouplingEfficiency(E_eff,sys_props,fiu_fiber_props,wvls,dx,N_iris)

figure;
plot(wvls*1e9,etas*100,'-o');
xlabel('Wavelength (nm)');
ylabel('Coupling efficiency (%)');
ylim([0 100])


%% Export data

% Write irradiance to a fits file
fitswrite(I_fp,['outputs/',outfile_label,'.fits']);

