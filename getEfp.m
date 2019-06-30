function [E,I] = getEfp(E_pup,E_support,sys_props,wvls,dx,dxi,N_pix)
%[E,I] = getEfp(E_pup,E_support,sys_props,wvls,dx,dxi,N_pix)
%   Function to generate an array containing the complex E-field and
%   irradiance distribution at the focal plane.
%
%   Inputs:
%       E_pup - Array (or cube) of complex field values at the pupil plane. 
%               dimensions are N_pad x N_pad x number of wavelengths. 
%       E_support - Array (or cube) of complex field values representing 
%                   an ideal wavefront.  
%       sys_props.f_lens2 - Focal length of the focussing lens [meters] 
%       wvls - vector of wavelengths [meters]
%       dx - Sample spacing [meters]
%       dxi - Sample spacing in the image plane, or pixel size [meters]
%       N_pix - Number of pixels in the final image 
%
%   Returns: 
%       E - Complex field (N_pix x N_pix x number of wavelengths)
%       I - Irradiance distribution. The no-mask PSF is peak normalized.

    % Get the input array size, and number of wavelengths
    [~,~,num_wavelengths] = size(E_pup);

    % E has dimensions N_pix x N_pix x number of wavelengths 
    E = zeros(N_pix,N_pix,numel(wvls)); % Empty cubes
    I = E;
    for wvl_index = 1:num_wavelengths % loop over wavelengths 
        
        wvl = wvls(wvl_index); % wavelength [meters] 
        
        E_pup_lam = E_pup(:,:,wvl_index); % get slice of cube
        E_support_lam = E_support(:,:,wvl_index);
        
        % Compute the E-field at the focal plane 
        E_fp0 = propcustom_mft_PtoF(E_support_lam,sys_props.f_lens2,wvl,dx,dxi,N_pix,dxi,N_pix);
        E_fp  = propcustom_mft_PtoF(E_pup_lam    ,sys_props.f_lens2,wvl,dx,dxi,N_pix,dxi,N_pix);

        normI = abs(E_fp0(N_pix/2+1,N_pix/2+1)).^2;
        E(:,:,wvl_index) = E_fp; % Add 2D E-field to the cube
        I(:,:,wvl_index) = abs(E_fp ).^2/normI;% Add 2D irradiance to the cube

    end
end

