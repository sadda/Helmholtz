function RieszGradient = ComputeRieszGradient(Gradient, TriInfo, matrices)
    freeIndices       = sum(TriInfo.phiProlongationVector,2) == 0;
    freeIndices6      = repmat(freeIndices, TriInfo.sizePhi, 1);
    GradientProlonged = ProlongPhi(Gradient, TriInfo) - TriInfo.phiProlongationVector(:);
    RieszGradientMod  = matrices.H1scal(freeIndices6,freeIndices6) \ GradientProlonged(freeIndices6);
    % Reduce the vector to half the size
    RieszGradient                        = zeros(TriInfo.sizePhi*TriInfo.npoint,1);
    RieszGradient(freeIndices6)          = RieszGradientMod;
    RieszGradient(~TriInfo.phiRowsFree6) = [];
end