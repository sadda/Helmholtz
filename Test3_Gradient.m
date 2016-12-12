% load('Results_Ref1_4_Helm/DataAll.mat');

iterMax = 1;
relErr  = zeros(4, iterMax);
sizePhi = TriInfo.sizePhi;

for i=1:iterMax
    [~,gradient]  = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,1);
    rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
    deltaPhi      = -reshape(rieszGradient,[],sizePhi);
    [err, DJ]     = Test3_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material);
    relErr(:,i)   = err(2:end)./DJ(2:end);
end

format long;
relErr


