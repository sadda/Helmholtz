function [J, G, J1, J21, J22, J3, G1, G21, G22, G3, u, Theta, dataEigen] = ComputeData2Mod(phi,TriInfo,Transformation,matrices,constants,material,computeG,dataEigen)
    
    regThetaL1 = 1;
    regPhiL1   = 1;
    regPhiL2   = 1;
    
    if nargin < 8 || isempty(dataEigen)
        dataEigen     = containers.Map('KeyType','double','ValueType','any');
        dataEigen(-1) = 0;
    end
    if nargin < 7
        computeG = 1;
    end
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    npoint   = TriInfo.npoint;
    idp      = TriInfo.idp;
    sizePhi  = TriInfo.sizePhi;
    
    constants.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    constants.epsilonR(2:3) = [];
    epsilonR                = constants.epsilonR;
    
    id1D = ~idp(1:npoint)==1;
    phi  = ProlongPhi(phi, TriInfo);
    
    %% Compute Theta and eigen
    
    phiSum = ComputePhiCutoffEigen(phi);
    phiSum = sum(phiSum.*repmat(epsilonR,npoint,1),2);
    % phi*Theta*v -> Theta
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
    T    = T(id1D,id1D);
    S    = matrices.GradSq(id1D,id1D);
    M1   = matrices.Mloc(1:npoint,1:npoint);
    M    = matrices.Mloc(id1D,id1D);
    
    if isKey(dataEigen, dataEigen(-1))
        data        = dataEigen(dataEigen(-1));
        yEigen      = data{1};
        lambdaEigen = data{2};
        phiEigen    = data{3};
        constM      = 1+1e-3;
        constL      = max(epsilonR);
        constCP     = 2*pi^2/((max(TriInfo.x)-min(TriInfo.x))*(max(TriInfo.y)-min(TriInfo.y)));
        constEigen  = constL*(constCP^2+1)*(2*constM*constCP^2+1)/(constCP^2);
        phiDiff     = phi - phiEigen;
        shift2      = -lambdaEigen + constEigen*sqrt(phiDiff(:)'*matrices.Mloc*phiDiff(:));
        shift2      = shift2 + 1e-5;
    else
        yEigen      = ones(sum(id1D),1);
        shift2      = Inf;
    end
    
    shift1  = max(epsilonR);
    shift   = min(shift1, shift2);
    maxIter = 2000;
    
    try
        R = chol(S-T+shift*M);
    catch
        R = chol(S-T+max(shift1,shift2)*M);
        warning('Shift was not successful. Switching for higher shift for this iteration.');
    end
    
    for i=1:maxIter
        x      = yEigen / norm(yEigen);
        Mx     = M*x;
        yEigen = R \ (R' \ Mx);
        eigen  = 1/(yEigen'*x);
        res    = (S-T+shift*M)*x-eigen*Mx;
        if norm(res) <= 1e-12
            break
        end
    end
    
    if i == maxIter
        if norm(res) <= 1e-6
            if norm(res) >= 1e-10
                warning('Maximum number of iterations for eigenvalue computation exceeded. The tolerance is still pretty low. This may happen at the beginning of projected gradients.');
            end
        else
            error('Maximum number of iterations for eigenvalue computation exceeded');
        end
    end
    Theta       = zeros(npoint,1);
    Theta(id1D) = x;
    eigen       = eigen - shift;
    Theta       = Theta / sqrt(Theta'*matrices.Mloc(1:npoint,1:npoint)*Theta);
    if mean(Theta) < 0
        Theta = -Theta;
    end
    
    
    
    
    
    if SymmetryError(Theta, TriInfo) >= 1e-8 || SymmetryError(Theta, TriInfo) >= 1e-8
        error('The symmetrization procedure failed');
    end
    yEigen = SymmetryCompute(yEigen, TriInfo);
    Theta  = SymmetryCompute(Theta, TriInfo);
    
    
    
    
    dataEigen(dataEigen(-1)) = {yEigen, eigen, phi};
    
    %% Compute objective
    
    J1  = 0.5*(Theta-TriInfo.phiGe)'*M1*(Theta-TriInfo.phiGe);
    J21 = ones(npoint,1)'*M1*Theta;
    J22 = ones(npoint,1)'*M1*phi(:,1);
    J3  = 0.5*phi(:,1)'*M1*phi(:,1);
    
    J   = J1 + regThetaL1*J21 + regPhiL1*J22 + regPhiL2*J3;
    u   = [];
    
    if nargout > 1 && computeG
        %% Compute q
        
        B11 = S-T-eigen*M;
        B12 = -M*Theta(id1D);
        B21 = (-M*Theta(id1D))';
        B22 = 0;
        
        B       = [B11 B12; B21 B22];
        bAdjQ1  = -M1*(Theta-TriInfo.phiGe);
        bAdjQ1  = [bAdjQ1(id1D); 0];
        bAdjQ21 = -M1*ones(npoint,1);
        bAdjQ21 = [bAdjQ21(id1D); 0];
        
        q1        = zeros(npoint, 1);
        q21       = zeros(npoint, 1);
        qr1       = B\bAdjQ1;
        qr21      = B\bAdjQ21;
        q1(id1D)  = qr1(1:end-1);
        q21(id1D) = qr21(1:end-1);
        
        %% Compute gradient
        
        % v*Theta*q -> q
        T   = zeros(nelement,1,1,9);
        for k=1:nelement
            % Cutoffs of g' are not needed here because we compute gradient only at feasible points, where no cutoffs take place
            edet   = Transformation{k,1};
            
            ThetaAux   = Theta(e2p(k,:));
            
            addition1  = 6*ThetaAux(1)+2*ThetaAux(2)+2*ThetaAux(3);
            addition5  = 2*ThetaAux(1)+6*ThetaAux(2)+2*ThetaAux(3);
            addition9  = 2*ThetaAux(1)+2*ThetaAux(2)+6*ThetaAux(3);
            addition24 = 2*ThetaAux(1)+2*ThetaAux(2)+1*ThetaAux(3);
            addition37 = 2*ThetaAux(1)+1*ThetaAux(2)+2*ThetaAux(3);
            addition68 = 1*ThetaAux(1)+2*ThetaAux(2)+2*ThetaAux(3);
            addition   = [addition1 addition24 addition37; addition24 addition5 addition68; addition37 addition68 addition9];
            
            for i=1:sizePhi
                T(k,i,1,:) = 1/120*epsilonR(i)*edet*addition(:);
            end
        end
        T  = sparse(TriInfo.ii3(:), TriInfo.jj3(:), T(:));
        
        G1            = -T*q1;
        G21           = -T*q21;
        G22           = zeros(size(G1));
        G22(1:npoint) = M1*ones(npoint,1);
        G3            = zeros(size(G1));
        G3(1:npoint)  = M1*phi(:,1);
        
        G1  = ShortenPhi(G1, TriInfo);
        G21 = ShortenPhi(G21, TriInfo);
        G22 = ShortenPhi(G22, TriInfo);
        G3  = ShortenPhi(G3, TriInfo);
                
        G  = G1 + regThetaL1*G21 + regPhiL1*G22 + regPhiL2*G3;
    else
        G  = [];
        G1 = [];
        G2 = [];
        G3 = [];
    end
end





function phi = ComputePhiCutoffEigen(phi)
    
    delta = 1e-3;
    
    set1  = phi<0;
    set2  = phi>1;
    phi(set1) = delta*atan(phi(set1)/delta);
    phi(set2) = 1 + delta*atan((phi(set2)-1)/delta);
end
