function [constants, material] = ObtainData(epsilon, alpha, wavelength, options)
    
    s         = 1/6;
    l         = 1e-8;
    constants = struct('epsilon', epsilon,'alpha',alpha,'s',s,'l',l);
    
    voidC         = epsilon^2;
    lSiO          =  16.071; mSiO =  20.798;
    lGe           =  44.279; mGe  =  27.249;
    lSiN          = 110.369; mSiN =  57.813;
    lvoid         = voidC*lSiO; mvoid=voidC*mSiO;
    lambda        = [lGe lSiN lSiO lvoid];
    mu            = [mGe mSiN mSiO mvoid];
    epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    epsilonR(2:3) = [];
    epsilonR      = (2*pi/wavelength)^2*epsilonR;
    if options.originalEps0
        eps0          = 2.5e-3;
        sigma0        = -2.5;
    else
        eps0          = 2.6e-3;
        sigma0        = -3.8;
    end
    material      = struct('lambda', lambda,'mu',mu,'eps0',eps0,'sigma0',sigma0,'epsilonR',epsilonR);
end