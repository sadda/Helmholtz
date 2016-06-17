iterMax = 5;
relErr = zeros(2, 3, iterMax);
sizePhi = TriInfo.sizePhi;
% Test random direction
for i=1:iterMax
    phi = rand(sum(TriInfo.phiRowsFree),sizePhi);
    phi = phi ./ repmat(sum(phi,2), 1, sizePhi);
    deltaPhi = rand(size(phi));
    [err, DJ] = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material);
    relErr(1,:,i) = err(2:end)./DJ(2:end);
end
squeeze(relErr(1,:,:))
% % Test direction of gradient
% for i=1:iterMax
%     phi = rand(size(phi));
%     phi = phi ./ repmat(sum(phi,2), 1, sizePhi);
%     [u,p]         = elasticity_adjoint(phi,TriInfo,Transformation,matrices,constants,material);
%     Gradient      = gradientJ(phi,u,p,TriInfo,Transformation,matrices,constants,material);
%     RieszGradient = ComputeRieszGradient(Gradient, TriInfo, matrices);
%     deltaPhi = -reshape(RieszGradient,[],sizePhi);
%     [err, DJ] = Test2_GradientPart(phi, deltaPhi, TriInfo, Transformation, matrices, constants, material);
%     relErr(2,:,i) = err(2:end)./DJ(2:end);
% end
% squeeze(relErr(1,:,:))
% squeeze(relErr(2,:,:))


