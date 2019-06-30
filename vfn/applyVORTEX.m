function Eout = applyVORTEX(Ein,vfn_props)
%Eout = applyVORTEX(Ein,vfn_props,coords)
%   Applies a vortex mask to the field 
%
%   Inputs:
%       Ein - Array (or cube) of complex field values. dimensions are 
%             N_pad x N_pad x number of wavelengths. 
%       vfn_props - struct with vfn properties. Example: 
%           vfn_props.charge = 2; 
%
%   Outputs:
%       Eout - Array (or cube) of complex field values. 
%              dimensions are N_pad x N_pad x number of wavelengths. 

    % Get the input array size, and number of wavelengths
    [Ngrid,~,num_wavelengths] = size(Ein);
    
    % Creates arrays with coordinates 
    coords = generateCoordinates(Ngrid);
    
    mask = generateVortexMask( vfn_props.charge, coords, [0 0] );
    
	% Eout has dimensions N_pad x N_pad x number of wavelengths 
    Eout = zeros(Ngrid,Ngrid,num_wavelengths); % Empty cube 
    for wvl_index = 1:num_wavelengths % loop over wavelengths 
        Ein_lam = Ein(:,:,wvl_index); % get slice of cube
        E_lam = Ein_lam.*mask; % propagate
        Eout(:,:,wvl_index) = E_lam; % insert slice of cube into output 
    end
end

