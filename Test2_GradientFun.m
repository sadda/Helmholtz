function Test2_GradientFun(phi,TriInfo,Transformation,matrices,constants,material,options)
    
    options.cutoffs         = 0;
    options.computeG        = 1;
    [J,G,J1,J2,J3,G1,G2,G3] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options);
    rieszGradient           = ComputeRieszGradient(G, TriInfo, matrices);
    deltaPhi                = -reshape(rieszGradient,[],TriInfo.sizePhi);
    
    err  = zeros(8,4);
    JAll = [J J1 J2 J3];
    GAll = [G G1 G2 G3];
    DJ   = zeros(1,4);
    for j=1:4
        DJ(j) = GAll(:,j)'*deltaPhi(:);
    end
    
    for k=1:10
        tau                = 10^(-k);
        phiNew             = phi+tau*deltaPhi;
        options.computeG   = 0;
        [JS,~,JS1,JS2,JS3] = ComputeData(phiNew,TriInfo,Transformation,matrices,constants,material,options);
        JDiff              = [JS JS1 JS2 JS3] - JAll;
        err(k,:)         = abs(JDiff/tau - DJ);
    end
    err = min(err);
    
    err ./ DJ
end



