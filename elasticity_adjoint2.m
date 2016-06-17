function [ThetaFull, lambda, adjoint1Full, adjoint2] = elasticity_adjoint2(phi, u, TriInfo, Transformation, matrices)
    
    npoint                  = TriInfo.npoint;
    nelement                = TriInfo.nelement;
    e2p                     = TriInfo.e2p;
    id                      = ~TriInfo.idp(1:TriInfo.npoint)==1;
    sizePhi                 = TriInfo.sizePhi;
    constants.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    constants.epsilonR(2:3) = [];
    epsilonR                = constants.epsilonR;
    
    %% Compute Theta and lambda
    
    phiPr  = ProlongPhi(phi,TriInfo);
    phiSum = sum(phiPr.*repmat(epsilonR,npoint,1),2);
    % phi*Theta*v
    T = zeros(nelement,1,1,9);
    for k=1:nelement
        edet   = Transformation{k,1};
        
        phiAux = phiSum(e2p(k,:));
        
        addition1  = 6*phiAux(1)+2*phiAux(2)+2*phiAux(3);
        addition5  = 2*phiAux(1)+6*phiAux(2)+2*phiAux(3);
        addition9  = 2*phiAux(1)+2*phiAux(2)+6*phiAux(3);
        addition24 = 2*phiAux(1)+2*phiAux(2)+1*phiAux(3);
        addition37 = 2*phiAux(1)+1*phiAux(2)+2*phiAux(3);
        addition68 = 1*phiAux(1)+2*phiAux(2)+2*phiAux(3);
        addition   = [addition1 addition24 addition37; addition24 addition5 addition68; addition37 addition68 addition9];
        
        T(k,1,1,:) = 1/120*edet*addition(:);
    end
    
    ii1D = TriInfo.indicesIPhi(:,1,1,:);
    jj2D = TriInfo.indicesJPhi(:,1,1,:);
    T    = sparse(ii1D(:), jj2D(:), T(:));
    T    = T(id,id);
    S    = matrices.GradSq(id,id);
    M    = matrices.Mloc(id,id);
    
    A = M \ (S-T);
    % A = diag(1./sum(M,2)) * (S-T);
    
    
    shift           = max(epsilonR);
    [Theta, lambda] = eigs(A + shift*speye(size(A)), 1, 'sm');
    lambda          = lambda - shift;
    Theta           = Theta / sqrt(Theta'*M*Theta);
    if median(Theta) < 0
        Theta = -Theta;
    end
    ThetaFull       = zeros(npoint,1);
    ThetaFull(id)   = Theta;
    
    %% Compute objective
    
    J = zeros(nelement,1,1,9);
    for k=1:nelement
        edet   = Transformation{k,1};
        slocx  = Transformation{k,9};
        slocy  = Transformation{k,10};
        mloc   = Transformation{k,11};
        
        ux   = u(e2p(k,:));
        uy   = u(e2p(k,:)+npoint);
        
        dxux = ux'*slocx';
        dyuy = uy'*slocy';
        
        J(k,1,1,:) = 2*(dxux+dyuy)/edet*mloc(:);
    end
    ii1D = TriInfo.indicesIPhi(:,1,1,:);
    jj2D = TriInfo.indicesJPhi(:,1,1,:);
    J    = sparse(ii1D(:), jj2D(:), J(:));
    J    = J(id,id);
    
    %% Adjoint equation
    
    B11 = S-T-lambda*M;
    B12 = -M*Theta;
    B21 = (M*Theta)';
    B22 = 0;
    
    B   = [B11 B12; B21 B22];
    b   = [2*J*Theta; 0];
    
    adjoint          =  B\b;
    adjoint1         = adjoint(1:end-1);
    adjoint1Full     = zeros(npoint,1);
    adjoint1Full(id) = adjoint1;
    adjoint2         = adjoint(end);
end












