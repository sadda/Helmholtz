function [err, DJ] = Test3_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material)    
    
    [J,G,J1,J21,J22,J3,G1,G21,G22,G3] = ComputeData2Mod(phi,TriInfo,Transformation,matrices,constants,material,1);    
    err  = zeros(16,5);
    JAll = [J J1 J21 J22 J3];
    GAll = [G G1 G21 G22 G3];
    for i=1:5
        DJ(i) = GAll(:,i)'*deltaPhi(:);
    end
    for k=0:15
        tau                      = 10^(-k);
        phiNew                   = phi+tau*deltaPhi;
        [JS,~,JS1,JS21,JS22,JS3] = ComputeData2Mod(phiNew,TriInfo,Transformation,matrices,constants,material,0);
        JDiff                    = [JS JS1 JS21 JS22 JS3] - JAll;
        err(k+1,:)               = abs(JDiff/tau - DJ);
    end
    err = min(err);
end