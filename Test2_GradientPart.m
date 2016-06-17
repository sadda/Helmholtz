function [err, DJ] = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material)    
    
    
    dummy = zeros(TriInfo.npoint,1);
    
    u                  = elasticity_adjoint(phi,dummy,TriInfo,Transformation,matrices,constants,material);
    [Theta,~,adjoint1] = elasticity_adjoint2(phi, u, TriInfo, Transformation, matrices);
    [~,p]              = elasticity_adjoint(phi,Theta,TriInfo,Transformation,matrices,constants,material);
    [J, J1, J2, J3]    = functionValue(phi,u,Theta,TriInfo,Transformation,matrices,constants);
    [G, G1, G2, G3]    = gradientJ(phi,u,p,Theta,adjoint1,TriInfo,Transformation,matrices,constants,material);
    %% compute approximation of gradient and error
    err=zeros(16,4);
    JAll = [J J1 J2 J3];
    GAll = [G G1 G2 G3];
    for i=1:4
        DJ(i) = GAll(:,i)'*deltaPhi(:);
    end
    for k=0:15
        tau                 = 10^(-k);
        phiNew              = phi+tau*deltaPhi;
        
        % TODO toto je uplne na hovno. protoze p zavisi na Theta.
        
        uNew                = elasticity_adjoint(phiNew,dummy,TriInfo,Transformation,matrices,constants,material);
        ThetaNew            = elasticity_adjoint2(phiNew, uNew, TriInfo, Transformation, matrices);
        
%         uNew = u;
%         ThetaNew = Theta;
        
        
        [JS, JS1, JS2, JS3] = functionValue(phiNew,uNew,ThetaNew,TriInfo,Transformation,matrices,constants);
        JDiff               = [JS JS1 JS2 JS3] - JAll;
        err(k+1,:)          = abs(JDiff/tau - DJ);
    end
    err = min(err);
end