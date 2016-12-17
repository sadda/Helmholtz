function err = SymmetryError(f, TriInfo, isPhi)

    fSym = SymmetryCompute(f, TriInfo, isPhi);
    err  = max(abs(fSym(:)-f(:)));
end