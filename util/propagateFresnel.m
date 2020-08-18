function Eout = propagateFresnel(Ein,z,wvls,dx)
%Eout = propagateFresnel(Ein,z,wvls,dx)
%   Propagates the E-field described by Ein a distance z using the Fresnel
%   approximation and angular spectrum method. 
%
%   Inputs:
%       Ein - Array (or cube) of complex field values. dimensions are 
%             N_pad x N_pad x number of wavelengths. 
%       z - Propagation distance [meters]
%       wvls - Wavelength(s) [meters] 
%       dx - Sample spacing [meters]
%
%   Outputs:
%       Eout - Array (or cube) of complex field values after propagation. 
%              dimensions are N_pad x N_pad x number of wavelengths. 


    if(z==0)
        Eout = Ein; 
    else
        % Get the input array size, and number of wavelengths
        [Ngrid,~,num_wavelengths] = size(Ein);

        % Creates arrays with coordinates 
        coords = generateCoordinates(Ngrid);

        % Fresnel propagation transfer function 
        H = @(lam) fftshift(exp(-1i*pi*lam*z*(coords.RHO/Ngrid/dx).^2));

        % Eout has dimensions N_pad x N_pad x number of wavelengths 
        Eout = zeros(Ngrid,Ngrid,num_wavelengths); % Empty cube 
        for wvl_index = 1:num_wavelengths % loop over wavelengths 
            wvl = wvls(wvl_index); % wavelength [meters] 
            Ein_lam = Ein(:,:,wvl_index); % get slice of cube
            E_lam = ifft2(fft2(Ein_lam).*H(wvl)); % propagate
            Eout(:,:,wvl_index) = E_lam; % insert slice of cube into output 
        end
        
    end
end

