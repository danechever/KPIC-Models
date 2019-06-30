function Eout = applyAPOD(Ein,apod_props,dx)
%Eout = applyAPOD(Ein,apod_props,wvls,dx)
%   Applies the apodizer to the complex field, Ein. The apodizer can be
%   loaded from a fits file or recalculated from coefficients. When using a
%   fits file, it is the users responsibility to ensure the samping matches
%   Ein. When using the cofficients, the apodization profile is a
%   polynomial and not error diffused. 
%
%   Inputs:
%       Ein - Array (or cube) of complex field values. dimensions are 
%             N_pad x N_pad x number of wavelengths. 
%       apod_props - struct with apodizer properties. Example: 
%           apod_props.D = sys_props.D_I;% diameter of the iris [meters]
%           apod_props.type = 'fits' or 'coeffs'; % Load from fits file. 
%           apod_props.whichdesign = 'KAPOD-01'; % Needed to choose fits file
%       dx - Sample spacing [meters]
%
%   Outputs:
%       Eout - Array (or cube) of complex field values. 
%              dimensions are N_pad x N_pad x number of wavelengths. 

    show_plots = true;
    
    % Get the input array size, and number of wavelengths
    [N_grid,~,num_wavelengths] = size(Ein);
    
    % Creates arrays with coordinates 
    coords = generateCoordinates(N_grid);
    
    if(strcmpi(apod_props.type,'fits'))
        try
            APOD0 = fitsread(['apod/',apod_props.whichdesign,'.fits']);
        catch
            gunzip(['apod/',apod_props.whichdesign,'.fits.gz']);
            APOD0 = fitsread(['apod/',apod_props.whichdesign,'.fits']);
        end
        APOD = padOrCropEven(APOD0,N_grid);
    else
        load('apod/coeff.mat','bestcoeffs','a');
        APOD = polyval(bestcoeffs,coords.RHO/(apod_props.D/dx/2));
        
        ROI = logical(makeKeckPupil( apod_props.D/dx, N_grid ));
        APOD = APOD - min(APOD(ROI));
        APOD = APOD/max(APOD(ROI));
        APOD(APOD<0) = 0;
        APOD(APOD>1) = 1;
    end
    
	% Eout has dimensions N_pad x N_pad x number of wavelengths 
    Eout = zeros(N_grid,N_grid,num_wavelengths); % Empty cube 
    for wvl_index = 1:num_wavelengths % loop over wavelengths 

        Ein_lam = Ein(:,:,wvl_index); % get slice of cube
        Eout_lam = Ein_lam.*APOD;
        Eout(:,:,wvl_index) = Eout_lam; % insert slice into cube
        
        if(show_plots)
            xvals = coords.xvals;yvals = coords.yvals;
            figure(104);imagesc(xvals*dx*1e3,yvals*dx*1e3,abs(Eout(:,:,wvl_index)));
                colorbar;axis image;set(gca,'ydir','normal');
                xlabel('x [mm]');ylabel('y [mm]');title('|E| after apodizer');
        end
        
    end
    
end

