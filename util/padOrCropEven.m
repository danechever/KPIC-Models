% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to pad or crop the even-sized input matrix to the desired size.
%
%--INPUTS:
% Ain: rectangular or square, even-sized matrix to be padded or cropped
% Ndes: Number of points desired across the array. Even number of points desired.
%
%--OUTPUTS:
% Aout: square, even-sized, padded or cropped matrix
%
% Modified by A.J. Riggs on 2018-03-06 to not rely on padarray.m. Dropping
% support for rectangular matrices right now for simplicity
% Modified by A.J. Riggs on 2017-10-27 to use varargin to allow the
%  extrapolation value to be specified.
% Can now pad or crop starting with rectangular matrices going to a square
%  matrix. -- A.J. Riggs, July 26, 2016

function Aout = padOrCropEven(Ain,Ndes,varargin)

    % Set default values of input parameters
    extrapval = 0; %--Value to use for extrapolated points

    %--Enable different extrapolation values by using varargin
    icav = 0; % index in cell array varargin
    while icav < size(varargin, 2)
        icav = icav + 1;
        switch lower(varargin{icav})
          case {'extrapval'}
            icav = icav + 1;
            extrapval   = varargin{icav};  %--Value to use for extrapolated points
          otherwise
            error('padOrCropEven: Unknown keyword: %s\n', ...
              varargin{icav});
        end
    end

    Nx0 = size(Ain,2);
    Ny0 = size(Ain,1);
    
	if(isa(Ain,'gpuArray'))
        Aout = gpuArray.ones(Ndes);
    end

    if( mod(Nx0,2)~=0 || mod(Ny0,2)~=0 )
        error('padOrCropEven: Input was not an even-sized array')
    elseif( Nx0~=Ny0 )
        error('padOrCropEven: Input was not a square array')
    elseif(length(Ndes)~=1) 
        error('padOrCropEven.m: Wrong number of dimensions specified for output')       
    else %--Pad or crop:
        if(min(Nx0,Ny0)>Ndes) %--Crop
            Aout = Ain( (Ny0-Ndes)/2+1:(Ny0+Ndes)/2, (Nx0-Ndes)/2+1:(Nx0+Ndes)/2 );
        elseif(max(Nx0,Ny0)<Ndes) %--Pad
            if(~isa(Ain,'gpuArray'))
                Aout = ones(Ndes); %--Initialize
            end
            Aout = extrapval*Aout;
            Aout( (Ndes-Ny0)/2+1:(Ndes+Ny0)/2, (Ndes-Nx0)/2+1:(Ndes+Nx0)/2 ) = Ain;
        else %--Leave as size
            Aout = Ain;
        end
    end
end %--END OF FUNCTION