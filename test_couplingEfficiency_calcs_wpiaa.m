%% test_couplingEfficiency_calcs.m 
% A script to test/demo the coupling efficiency calculations 

clear ;

%% Inputs 

% Sampling parameters
Nbeam = 500;
N_lambdaFnum = 4;

wvls = (2:0.1:4.0)*1e-6;% wavelengths [meters]
wvl_design = 3.7e-6; % design wavelength [meters]

%%- Fiber properties (properties store in fiber_props struct)

% % Parameters for Thorlabs SM2000
% % link: https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=949&pn=SM2000#4668
% fiber_props.core_rad = 11e-6/2;% Core radius [meters]
% fiber_props.n_core = 1.4436;% core index
% fiber_props.n_clad = 1.4381;% cladding index

% Parameters for the ZBLAN fiber 
fiber_props.core_rad = 6.5e-6;% Core radius [meters]
fiber_props.NA = 0.175;% Measured to be 0.175+-.01

% fiber_props.type = 'gaussian';
fiber_props.type = 'bessel'; % This one is more accurate 

% sys_props.pupil_shape = 'circ';
% sys_props.pupil_shape = 'keck';
sys_props.pupil_shape = 'kecklab';

sys_props.D_I = 12.3e-3; % diameter of the input beam [meters] 

% Set the focal length 
% use the optimal focal length at the design wavelength 
% sys_props.f_fiber = getMFD(fiber_props,wvl_design)*sys_props.D_I/wvl_design/1.4;
sys_props.f_fiber = 36.6e-3;% alternatively, use the actual focal length

%%- PIAA properties 
piaa_props.filename = 'piaa/KPIAA-01.csv';% file containing sag profiles 
piaa_props.L = 50e-3; % Distance between PIAA surfaces [meters]
piaa_props.material = 'CaF2'; % piaa lens material 
piaa_props.Dstop = 20e-3; % Diameter of the output beam [meters]

apply_piaa = true; 

%% Make the pupil mask 

%%- Make a 2D array containing the iris 
switch lower(sys_props.pupil_shape)
    case 'circ'
        PUPIL = makeCircularPupil(Nbeam/2, 2*Nbeam );
    case 'keck'
        [PUPIL,Nbeam] = makeKeckPupil( Nbeam, 2*Nbeam );
    case 'kecklab'
        [PUPIL,Nbeam] = makeKeckLabPupil( Nbeam, 2*Nbeam );
    otherwise
        error('sys_prop.pupil_shape not recognized.');
end
Narr = 2^nextpow2(Nbeam*N_lambdaFnum);
PUPIL = padOrCropEven(PUPIL,Narr);
dx = sys_props.D_I/Nbeam; % sample spacing [meters]

% E has dimensions Narr x Narr x number of wavelengths 
E = zeros(Narr,Narr,numel(wvls)); % Empty cube 
for wvl_index = 1:numel(wvls) % loop over wavelengths 
    E(:,:,wvl_index) = PUPIL; % Add 2D E-field to the cube
end

if(apply_piaa)
    addpath('piaa')
    E = applyPIAA(E,piaa_props,wvls,dx);
end
%% 

etas = getCouplingEfficiency(E,sys_props,fiber_props,wvls,dx,Nbeam);

figure;
plot(wvls*1e9,etas*100,'-o');
xlabel('Wavelength (nm)');
ylabel('Coupling efficiency (%)');
ylim([0 100])

