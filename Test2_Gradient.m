iterMax = 1;
relErr  = zeros(2, 3, iterMax);
sizePhi = TriInfo.sizePhi;
% Test random direction
for i=1:iterMax
    phi           = rand(sum(TriInfo.phiRowsFree),sizePhi);
    phi           = phi ./ repmat(sum(phi,2), 1, sizePhi);
    deltaPhi      = rand(size(phi));
    [err, DJ]     = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material);
    relErr(1,:,i) = err(2:end)./DJ(2:end);
end
squeeze(relErr(1,:,:))
% Test direction of gradient
for i=1:iterMax
    phi           = rand(size(phi));
    phi           = phi ./ repmat(sum(phi,2), 1, sizePhi);
    [~,gradient]  = ComputeData(phi,TriInfo,Transformation,matrices,constants,material);
    rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
    deltaPhi      = -reshape(rieszGradient,[],sizePhi);
    [err, DJ]     = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material);
    relErr(2,:,i) = err(2:end)./DJ(2:end);
end
squeeze(relErr(1,:,:))
squeeze(relErr(2,:,:))


