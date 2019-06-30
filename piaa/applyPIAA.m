function Eout = applyPIAA(Ein,piaa_props,wvls,dx)
%Eout = applyPIAA(Ein,piaa,wvls,dx)
%   Computes the effect of the piaa lenses on the complex field using the
%   following steps:
%       1. apply the 2D phase shift due to the first lens.
%       2. propagate to the second.
%       3. apply the 2D phase shift due to the second lens.
%       4. propagate back to the original surface. 
%       5. Return the effective field at the first surface 
%
%   Inputs:
%       Ein - Array (or cube) of complex field values. dimensions are 
%             N_pad x N_pad x number of wavelengths. 
%       piaa_props - struct with piaa properties. Example: 
%           piaa_props.filename = 'piaa/KPIAA-01.csv';% file containing sag profiles 
%           piaa_props.L = 50e-3; % Distance between PIAA surfaces [meters]
%           piaa_props.material = 'CaF2'; % piaa lens material 
%           piaa_props.Dstop = 20e-3; % Diameter of the output beam [meters]
%       wvls - Wavelength(s) [meters] 
%       dx - Sample spacing [meters]
%
%   Outputs:
%       Eout - Array (or cube) of complex field values. 
%              dimensions are N_pad x N_pad x number of wavelengths. 

    show_plots = true;

    % Get piaa profiles from file 
    M = csvread(piaa_props.filename,1,0);
    piaaRvals = M(:,1)*1e-3;% meters
    piaaZlens1 = M(:,2)*1e-3;% meters
    piaaZlens2 = M(:,5)*1e-3;% meters
    L = piaa_props.L; % Distance between PIAA surfaces [meters]
    
    % Get the input array size, and number of wavelengths
    [N_grid,~,num_wavelengths] = size(Ein);
    
    % Creates arrays with coordinates 
    coords = generateCoordinates(N_grid);
    
    % Get 2D arrays containing PIAA shapes 
    sagLens1 = interp1(piaaRvals,piaaZlens1,coords.RHO*dx,'linear','extrap');
    sagLens2 = interp1(piaaRvals,piaaZlens2,coords.RHO*dx,'linear','extrap');

    if(show_plots)
        xvals = coords.xvals;yvals = coords.yvals;
        figure(101);imagesc(xvals*dx*1e3,yvals*dx*1e3,sagLens1);
            colorbar;axis image;set(gca,'ydir','normal');
            xlabel('x [mm]');ylabel('y [mm]');title('Sag 1 [mm]');
        figure(102);imagesc(xvals*dx*1e3,yvals*dx*1e3,sagLens2);
            colorbar;axis image;set(gca,'ydir','normal');
            xlabel('x [mm]');ylabel('y [mm]');title('Sag 2 [mm]');
    end
    
    % Aperture stop applied at second lens
    STOP = makeCircularPupil(piaa_props.Dstop/dx/2, N_grid );
    
	% Eout has dimensions N_pad x N_pad x number of wavelengths 
    Eout = zeros(N_grid,N_grid,num_wavelengths); % Empty cube 
    for wvl_index = 1:num_wavelengths % loop over wavelengths 
        
        wvl = wvls(wvl_index); % wavelength [meters] 
        
        % Define the refractive index of each lens 
        n1 = getRefractiveIndex(piaa_props.material,wvl*1e6);
        n2 = getRefractiveIndex(piaa_props.material,wvl*1e6);
        
        phzLENS1 = 2*pi/wvl*(n1-1)*sagLens1;
        phzLENS2 = 2*pi/wvl*(1-n2)*sagLens2;
        
        Ein_lam = Ein(:,:,wvl_index); % get slice of cube

        E1 = Ein_lam.*exp(1i*phzLENS1);
        E2 = propagateFresnel(E1,L,wvls,dx);
        E3 = E2.*exp(1i*phzLENS2).*STOP;
        E4 = propagateFresnel(E3,-L,wvls,dx);

        Eout(:,:,wvl_index) = E4; % insert slice into cube
        
        if(show_plots)
            figure(103);imagesc(xvals*dx*1e3,yvals*dx*1e3,abs(E2));
                colorbar;axis image;set(gca,'ydir','normal');
                xlabel('x [mm]');ylabel('y [mm]');title('|E| at PIAA lens 2');
        end
    
        
    end

end

