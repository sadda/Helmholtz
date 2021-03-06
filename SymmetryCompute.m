function [fSym, symError] = SymmetryCompute(f, TriInfo, isPhi, writeError, stopError)
    % Symmetrizes phi or Theta
    
    %% Prolongs the function onto the whole mesh
    
    fOrig = f;
    if isPhi && size(f,1) > TriInfo.npoint
        f = reshape(f, [], TriInfo.sizePhi);
    end
    fReshaped = f;
    % Reduced phi
    if isPhi && size(f,1) == TriInfo.npointRed
        fAux = zeros(TriInfo.npoint, size(f,2));
        for i=1:size(f,2)
            fAux(:,i) = TriInfo.phiProlongationMatrix*f(:,i);
        end
        f = fAux;
    end
    % Theta
    if ~isPhi && size(f,1) == sum(~TriInfo.idp)
        fAux = zeros(TriInfo.npoint, size(f,2));
        fAux(~TriInfo.idp,:) = f;
        f = fAux;
    end
    
    %% Finds the node symmetric along the x axis and computes the mean valus of these two
    
    fSym = f;
    for i=1:size(f,1)
        if TriInfo.x(i) < 0
            ind = find(TriInfo.x == -TriInfo.x(i) & TriInfo.y == TriInfo.y(i));
            if length(ind) == 1
                fSym(i,:)   = 0.5*(f(i,:)+f(ind,:));
                fSym(ind,:) = 0.5*(f(i,:)+f(ind,:));
            else
                error('Something wrong with the mesh');
            end
        end
    end
    
    %% Transforms the symmetrized functions into the desired shape
    
    if isPhi && size(fReshaped,1) == TriInfo.npointRed
        fSym = fSym(TriInfo.phiRowsFree,:);
    end
    if ~isPhi && size(fReshaped,1) == sum(~TriInfo.idp)
        fSym = fSym(~TriInfo.idp,:);
    end
    if size(fOrig,1) > TriInfo.npoint
        fSym = fSym(:);
    end
    
    %% Computes the symmetrization error and possibly write a warning if it too large
    
    symError = norm(fSym - fOrig);
    if nargin >= 4 && writeError
        fprintf('\nThe symmetrization error was %1.3e.\n', symError) ;
    end
    if nargin >= 5 && ~isempty(stopError) && symError > stopError
        warning('The symmetrization error was too big: %1.3e.\n', symError);
    end
end



