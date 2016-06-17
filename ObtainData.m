function [constants, material] = ObtainData(epsilon, alpha)
    s         = 1/6;
    l         = 1e-8;
    constants = struct('epsilon', epsilon,'alpha',alpha,'s',s,'l',l);
    
    voidC  = epsilon^2;
    lSiO   =  16.071; mSiO =  20.798; % SiO2
    % lSip   =  64.63 ; mSip =  33.855; % Si_p
    % lSin   =  64.63 ; mSin =  33.855; % Si_n
    lGe    =  44.279; mGe  =  27.249; % Ge
    lSiN   = 110.369; mSiN =  57.813; % SiN
    lvoid  = voidC*lSiO; mvoid=voidC*mSiO;    % void
    lambda = [lGe lSiN lSiO lvoid];
    mu     = [mGe mSiN mSiO mvoid];
    
    eps0   = 2.5e-3;
    sigma0 = -2.5;
    material = struct('lambda', lambda,'mu',mu,'eps0',eps0,'sigma0',sigma0);
end