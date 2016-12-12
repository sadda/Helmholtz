function err = SymmetryError(f, TriInfo)
    
    fSym = SymmetryCompute(f, TriInfo);
    err  = max(abs(fSym(:)-f(:)));
end