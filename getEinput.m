function [E,coords,w_I] = getEinput(fiber_props,sys_props,wvls,N_iris,N_pad)
%[E,coords,w_I] = getEinput(fiber_props,sys_props,wvls,N_iris,N_pad)
%   Function to generate an array containing the complex E-field in a
%   system where a single-mode fiber is the source and the optical system
%   is descibed by a ray transfer matrix or 'ABCD' matrix (see e.g. 
%   https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis). 
%   The ray transfer matrix should describe the propagation from the fiber
%   tip to the plane of interest. The returned E-field is peak normalized. 
%
%   Inputs: 
%       - Fiber definition - 
%       fiber_props.core_rad - Core radius of the fiber [meters]
%       fiber_props.NA - Numerical aperture of the fiber (optional)
%       * If NA is omitted, the following optional parameters are used instead *
%       fiber_props.n_core - Refractive index of the fiber core
%       fiber_props.n_clad - Refractive index of the fiber cladding
%
%       - System definition - 
%       sys_props.TM - 2x2 ray transfer matrix (fiber -> plane of interest)
%       sys_props.D_I - Diameter of the iris [meters]
%
%       - Beam properties - 
%       wvls - vector of wavelengths [meters]
%       N_iris - Number of samples across the iris 
%       N_pad - Number of samples across the full array 
%
%   Returns: 
%       E - Complex field (N_pad x N_pad x number of wavelengths)
%       coords - coordinate system in units of samples (optional) 
%       w_I - Beam waist at the plane of interest (optional) 

    %%- Get beam size at the source fiber tip
    MFD = getMFD(fiber_props,wvls);% mode field diameters at the specified wavelengths 
    w_F = MFD/2; % Beam waist at the fiber [meters]
    z_R = pi*w_F.^2./wvls; % Rayleigh range [meters]
    
    %%- Unpack system inputs 
    TM = sys_props.TM;  % 2x2 transfer (or ABCD) matrix
    D_I = sys_props.D_I;% diameter of the iris [meters]

    %%- compute the complex beam parameters as a function of wavelength
    % see e.g. https://en.wikipedia.org/wiki/Gaussian_beam
    q_F = 1i*z_R;% q at fiber 
    q_I_inv = (TM(2)+TM(4)./q_F)./(TM(1)+TM(3)./q_F);% 1/q at iris
    w_I = abs(sqrt(wvls./(pi*imag(q_I_inv))));% beam waist at the iris

    %%- Make a 2D array containing the iris 
    switch lower(sys_props.pupil_shape)
        case 'circ'
            IRIS = makeCircularPupil(N_iris/2, 2*N_iris );
        case 'keck'
            [IRIS,N_iris] = makeKeckPupil( N_iris, 2*N_iris );
        case 'kecklab'
            [IRIS,N_iris] = makeKeckLabPupil( N_iris, 2*N_iris );
        otherwise
            error('sys_prop.pupil_shape not recognized. options: circ or keck');
    end
    IRIS = padOrCropEven(IRIS,N_pad);
    
    %%- Build the E-field cube at the plane of interest 
    coords = generateCoordinates(N_pad);% Creates arrays with coordinates 
    dx = D_I/N_iris; % sample spacing at iris [meters]
    
    % E has dimensions N_pad x N_pad x number of wavelengths 
    E = zeros(N_pad,N_pad,numel(wvls)); % Empty cube 
    for wvl_index = 1:numel(wvls) % loop over wavelengths 
        wvl = wvls(wvl_index); % wavelength [meters] 
        E_lam = exp(-1i*pi*(coords.RHO*dx).^2.*q_I_inv(wvl_index)/wvl); % E-field
        E(:,:,wvl_index) = E_lam.*IRIS; % Add 2D E-field to the cube
    end
    
end

