function [err, DJ] = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material)    
    
    [J,G,J1,J2,J3,G1,G2,G3] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material);    
    err  = zeros(16,4);
    JAll = [J J1 J2 J3];
    GAll = [G G1 G2 G3];
    for i=1:4
        DJ(i) = GAll(:,i)'*deltaPhi(:);
    end
    for k=0:15
        tau                = 10^(-k);
        phiNew             = phi+tau*deltaPhi;
        [JS,~,JS1,JS2,JS3] = ComputeData(phiNew,TriInfo,Transformation,matrices,constants,material,0);
        JDiff              = [JS JS1 JS2 JS3] - JAll;
        err(k+1,:)         = abs(JDiff/tau - DJ);
    end
    err = min(err);
end