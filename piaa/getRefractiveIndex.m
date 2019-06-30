function n = getRefractiveIndex(material,lambda)
%n = getRefractiveIndex(material,lambda)
%   Returns refractive index of 'material' at wavelength lambda.
%   
%   Inputs:
%       material - String indicating which material (e.g. 'CaF2')
%       lambda - Wavelength in microns 
%
%   Outputs:
%       n - (Real) refractive index

    if(material == 'CaF2')
        % Calcium Floride Sellmeier coefficients from Thorlabs 
        % Originally from H. H. Li., J. Phys. Chem. Ref. Data 9, 161-289 (1980)
        % https://aip.scitation.org/doi/pdf/10.1063/1.555616
        C0=0.33973;
        B1=0.69913;
        C1=0.09374^2;
        B2=0.11994;
        C2=21.18^2;
        B3=4.35181;
        C3=38.46^2;
    else
        error([material,' not yet implemented.']);
    end
    
    n = sqrt(1 + C0 + B1*lambda.^2./(lambda.^2-C1) + ...
                          B2*lambda.^2./(lambda.^2-C2) + ...
                          B3*lambda.^2./(lambda.^2-C3));
                      
end

