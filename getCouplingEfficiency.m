function eta = getCouplingEfficiency(E_pup,sys_props,fiber_props,wvls,dx,Nbeam)
%eta = getCouplingEfficiency(E_pup,E_support,sys_props,fiber_props,wvls,dx)
%   Function to generate coupling efficiency versus wavelength 
%
%   Inputs:
%       E_pup - Array (or cube) of complex field values at the pupil plane. 
%               dimensions are N_pad x N_pad x number of wavelengths. 
%       E_support - Array (or cube) of complex field values representing 
%                   an ideal wavefront.  
%       sys_props.f_fiber - Focal length of the focussing lens [meters] 
%       fiber_props - Fiber properties 
%       wvls - vector of wavelengths [meters]
%       dx - Sample spacing [meters]
%
%   Returns: 
%       eta - Coupling efficiency for each wavelength in wvls

    % Get the input array size, and number of wavelengths
    [~,~,num_wavelengths] = size(E_pup);

    dxi = fiber_props.core_rad/10;
    N = 200; % computes over 200x200 samples 
    coords = generateCoordinates(N);
    
    for wvl_index = 1:num_wavelengths % loop over wavelengths 
        
        wvl = wvls(wvl_index); % wavelength [meters] 
        
        E_pup_lam = E_pup(:,:,wvl_index); % get slice of cube
        E_pup_lam = padOrCropEven(E_pup_lam,2*round(Nbeam));
        
        totalPower = sum(abs(E_pup_lam(:).^2));
        
        % Compute the E-field at the focal plane 
        E_fp  = propcustom_mft_PtoF(E_pup_lam,sys_props.f_fiber,wvl,dx,dxi,N,dxi,N);
        
        % Generate the fiber mode 
        Fnum = sys_props.f_fiber/(dx*Nbeam);
        lambdaFnum_samps = wvl*Fnum/dxi;
        mode = generateSMFmode( fiber_props, wvl, Fnum, lambdaFnum_samps, coords);

        % Plots for debugging 
%         figure(11);imagesc(abs(E_fp).^2/max(abs(E_fp(:)).^2));
%         figure(12);imagesc(mode);
        figure(13);
        plot(abs(E_fp(N/2+1,:)).^2/max(abs(E_fp(:)).^2));hold on;
        plot(mode(N/2+1,:)/max(mode(:)));hold off;
        
        % Calculate the coupling efficiency 
        eta(wvl_index) = abs(sum(sum(mode.*E_fp))).^2/totalPower;

    end
end

