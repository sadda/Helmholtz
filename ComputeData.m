function [J, G, J1, J2, J3, G1, G2, G3, u, Theta, dataEigen] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,computeG,dataEigen)
    
    if nargin < 7
        computeG = 1;
    end
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    npoint   = TriInfo.npoint;
    idp      = TriInfo.idp;
    sizePhi  = TriInfo.sizePhi;
    
    lambda   = material.lambda;
    mu       = material.mu;
    eps0     = material.eps0;
    sigma0   = material.sigma0;
    
    alpha    = constants.alpha;
    epsilon  = constants.epsilon;
    
    constants.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    constants.epsilonR(2:3) = [];
    epsilonR                = constants.epsilonR;
    
    id             = ~[idp(1:npoint)==1;idp(1:npoint)==1];
    id1D           = ~[idp(1:npoint)==1];
    phi            = ProlongPhi(phi, TriInfo);
    
    %% Compute u
    slocx3a_aa = matrices.slocx3a_aa;
    slocxx_aa  = matrices.slocxx_aa;
    slocy3a_aa = matrices.slocy3a_aa;
    slocyy_aa  = matrices.slocyy_aa;
    slocxy_aa  = matrices.slocxy_aa;
    slocyx_aa  = matrices.slocyx_aa;
    
    AEla11_2_1  = zeros(nelement,2,sizePhi,9);
    AEla12      = zeros(nelement,2,2,9);
    AEla12_11_1 = zeros(nelement,1);
    AEla12_11_2 = zeros(nelement,1);
    AEla12_22_1 = zeros(nelement,1);
    AEla12_22_2 = zeros(nelement,1);
    AEla12_12_1 = zeros(nelement,1);
    AEla12_12_2 = zeros(nelement,1);
    
    % C(phi)e(u):e(v) -> u
    for i=1:TriInfo.sizePhi
        phiAux = phi(:,i);
        phiAux = ComputePhiCutoff(phiAux);
        sumPhi = sum(phiAux(e2p),2);
        AEla12_11_1 = AEla12_11_1 + (2/3*mu(i)+1/3*lambda(i))*sumPhi;
        AEla12_11_2 = AEla12_11_2 + 1/3*mu(i)*sumPhi;
        AEla12_22_1 = AEla12_22_1 + 1/3*mu(i)*sumPhi;
        AEla12_22_2 = AEla12_22_2 + (2/3*mu(i)+1/3*lambda(i))*sumPhi;
        AEla12_12_1 = AEla12_12_1 + 1/3*mu(i)*sumPhi;
        AEla12_12_2 = AEla12_12_2 + 1/3*lambda(i)*sumPhi;
    end
    AEla12(:,1,1,:) = repmat(AEla12_11_1,1,1,1,9).*slocxx_aa + repmat(AEla12_11_2,1,1,1,9).*slocyy_aa;
    AEla12(:,2,2,:) = repmat(AEla12_22_1,1,1,1,9).*slocxx_aa + repmat(AEla12_22_2,1,1,1,9).*slocyy_aa;
    AEla12(:,1,2,:) = repmat(AEla12_12_1,1,1,1,9).*slocyx_aa + repmat(AEla12_12_2,1,1,1,9).*slocxy_aa;
    AEla12(:,2,1,:) = repmat(AEla12_12_1,1,1,1,9).*slocxy_aa + repmat(AEla12_12_2,1,1,1,9).*slocyx_aa;
    % phi_i*phi1_1*tr(e(v)) -> phi, phi1
    phiAux = phi(:,1);
    phiAux = phiAux(e2p);
    for i = 1:sizePhi
        aux1x   = 1/6*eps0*(mu(i)+lambda(i))*slocx3a_aa.*reshape(phiAux(:,[1 1 1 2 2 2 3 3 3]), [], 1, 1, 9);
        aux1y   = 1/6*eps0*(mu(i)+lambda(i))*slocy3a_aa.*reshape(phiAux(:,[1 1 1 2 2 2 3 3 3]), [], 1, 1, 9);
        aux2x   = 1/6*eps0*(mu(i)+lambda(i))*slocx3a_aa.*repmat(sum(phiAux,2), 1, 1, 1, 9);
        aux2y   = 1/6*eps0*(mu(i)+lambda(i))*slocy3a_aa.*repmat(sum(phiAux,2), 1, 1, 1, 9);
        AEla11_2_1(:,1,i,:) = aux1x + aux2x;
        AEla11_2_1(:,2,i,:) = aux1y + aux2y;
    end
    
    bEla = zeros(nelement,2,3);
    for k=1:nelement
        clocx = Transformation{k,7};
        clocy = Transformation{k,8};
        
        phiMatrix0ele = [phi(e2p(k,1),:)',phi(e2p(k,2),:)',phi(e2p(k,3),:)'];
        
        sigmaXphi     = sigma0*phiMatrix0ele(2,:);
        bEla( k,1,:) = sigmaXphi*clocx ;
        bEla( k,2,:) = sigmaXphi*clocy ;
    end
    
    AEla       = sparse(TriInfo.indicesIEla(:), TriInfo.indicesJEla(:), AEla12(:));
    AEla11_2_1 = sparse(TriInfo.ii1(:), TriInfo.jj1(:), AEla11_2_1(:));
    bEla       = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bEla(:));
    bEla       = AEla11_2_1*phi(:) - bEla;
    
    u          = zeros(2*npoint,1);
    u(id)      = AEla(id,id) \ bEla(id);
    
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
    
    ii1D  = TriInfo.indicesIPhi(:,1,1,:);
    jj2D  = TriInfo.indicesJPhi(:,1,1,:);
    T     = sparse(ii1D(:), jj2D(:), T(:));
    T     = T(id1D,id1D);
    S     = matrices.GradSq(id1D,id1D);
    M     = matrices.Mloc(id1D,id1D);
    
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
        yEigen      = rand(sum(id1D),1);
        shift2      = Inf;
    end
    
    shift1  = max(epsilonR);
    shift   = min(shift1, shift2);
    
    maxIter = 200;
    
    R       = chol(S-T+shift*M);
    Theta   = zeros(npoint,1);
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
        warning('Maximum number of iterations for eigenvalue computation exceeded');
    end
    Theta(id1D) = x;
    eigen       = eigen - shift;
    Theta       = Theta / sqrt(Theta'*matrices.Mloc(1:npoint,1:npoint)*Theta);
    if median(Theta) < 0
        Theta = -Theta;
    end
        
    dataEigen(dataEigen(-1)) = {yEigen, eigen, phi};
    
    %     eigenOpt.v0 = yEigen;
    %     qwe = eigs(S-T+shift*M, M, 1, 'sm', eigenOpt);
    %     qwe = qwe - shift;
    %     norm(qwe-eigen)
    
    %% Compute objective
    
    % tr(e(u))theta^2 -> theta, theta
    AObj = zeros(nelement,1,1,9);
    for k=1:nelement
        edet   = Transformation{k,1};
        slocx  = Transformation{k,9};
        slocy  = Transformation{k,10};
        mloc   = Transformation{k,11};
        
        ux     = u(e2p(k,:));
        uy     = u(e2p(k,:)+npoint);
        
        dxux   = ux'*slocx';
        dyuy   = uy'*slocy';
        
        AObj(k,1,1,:) = 2*(dxux+dyuy)/edet*mloc(:);
    end
    ii1D = TriInfo.indicesIPhi(:,1,1,:);
    jj2D = TriInfo.indicesJPhi(:,1,1,:);
    AObj = sparse(ii1D(:), jj2D(:), AObj(:));
    
    J1   = -Theta'*AObj*Theta;
    J2   = 0.5*(matrices.GradSq*phi(:))'*phi(:);
    J3   = 0.5*(matrices.Id'*phi(:) - (matrices.Mloc*phi(:))'*phi(:));
    
    J    = J1 + alpha*epsilon*J2 + alpha/epsilon*J3;
    
    if nargout > 1 && computeG
        %% Compute p
        
        % tr(e(v))theta^2 -> v (vector)
        bAdjP = zeros(nelement,1,3);
        for k=1:nelement
            edet   = Transformation{k,1};
            slocx  = Transformation{k,9};
            slocy  = Transformation{k,10};
            mloc   = Transformation{k,11};
            
            ThetaAux = Theta(e2p(k,:));
            
            bAdjP(k,1,:) = 2*slocx'/edet*(ThetaAux'*mloc*ThetaAux);
            bAdjP(k,2,:) = 2*slocy'/edet*(ThetaAux'*mloc*ThetaAux);
        end
        bAdjP   = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bAdjP(:));
        
        p       = zeros(2*npoint,1);
        p(id)   = AEla(id,id)\bAdjP(id);
        
        %% Compute q
        
        B11 = S-T-eigen*M;
        B12 = -M*Theta(id1D);
        B21 = (M*Theta(id1D))';
        B22 = 0;
        
        B     = [B11 B12; B21 B22];
        bAdjQ = [2*AObj(id1D,:)*Theta; 0];
        
        q       = zeros(npoint, 1);
        qr      = B\bAdjQ;
        q(id1D) = qr(1:end-1);
        
        %% Compute gradient
        
        EdPhi  = zeros(nelement,sizePhi,3);
        for k=1:nelement
            slocx = Transformation{k,9};
            slocy = Transformation{k,10};
            
            pxloc = [p(e2p(k,1)); p(e2p(k,2)); p(e2p(k,3))];
            pyloc = [p(e2p(k,1)+npoint); p(e2p(k,2)+npoint); p(e2p(k,3)+npoint)];
            
            dxpx = pxloc'*slocx';
            dypy = pyloc'*slocy';
            
            if sizePhi == 6
                EdPhi(k,4,:) = (sigma0*2*dxpx+sigma0*2*dypy)*1/6*ones(3,1);
            else
                EdPhi(k,2,:) = (sigma0*2*dxpx+sigma0*2*dypy)*1/6*ones(3,1);
            end
        end
        
        A41_4      = zeros(nelement,sizePhi,sizePhi,9);
        A42        = zeros(nelement,sizePhi,2,9);
        edet_aa    = matrices.edet_aa;
        mloc_aa    = matrices.mloc_aa;
        slocx_aa   = matrices.slocx_aa;
        slocx3b_aa = matrices.slocx3b_aa;
        slocy_aa   = matrices.slocy_aa;
        slocy3b_aa = matrices.slocy3b_aa;
        
        px     = p(1:npoint);
        py     = p(npoint+1:end);
        px     = px(TriInfo.e2p);
        py     = py(TriInfo.e2p);
        % (phi_1*w_i+phi_i*w_1)*tr(e(p)) -> phi
        aux1 = 8*(lambda(1)+mu(1))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
        A41_4(:,1,1,:) = repmat(aux1,1,1,1,9).*mloc_aa;
        for m=2:sizePhi
            aux2 = 4*(lambda(m)+mu(m))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
            A41_4(:,1,m,:) = repmat(aux2,1,1,1,9).*mloc_aa;
            A41_4(:,m,1,:) = repmat(aux2,1,1,1,9).*mloc_aa;
        end
        EdPhi1    = sparse(TriInfo.indicesIPhi(:), TriInfo.indicesJPhi(:), -A41_4(:));
        % C(w)e(u):e(p) -> u
        for m=1:sizePhi
            pxBackup    = px;
            pyBackup    = py;
            cDerivative = ComputeCDerivative(phi(:,m), TriInfo);
            px          = px.*cDerivative;
            py          = py.*cDerivative;
            
            aux1_1Aux = reshape(px,[],1,1,3).*slocx_aa + reshape(py,[],1,1,3).*slocy_aa;
            aux1_1    = 2/3*lambda(m)*repmat(sum(aux1_1Aux,4)./edet_aa,1,1,1,9).*slocx3b_aa;
            aux1_2Aux = reshape(px,[],1,1,3).*slocx_aa + reshape(py,[],1,1,3).*slocy_aa;
            aux1_2    = 2/3*lambda(m)*repmat(sum(aux1_2Aux,4)./edet_aa,1,1,1,9).*slocy3b_aa;
            aux2_1Aux = reshape(px,[],1,1,3).*slocx_aa;
            aux2_1    = 4/3*mu(m)*repmat(sum(aux2_1Aux,4)./edet_aa,1,1,1,9).*slocx3b_aa;
            aux2_2Aux = reshape(py,[],1,1,3).*slocy_aa;
            aux2_2    = 4/3*mu(m)*repmat(sum(aux2_2Aux,4)./edet_aa,1,1,1,9).*slocy3b_aa;
            aux2_3Aux = reshape(py,[],1,1,3).*slocx_aa + reshape(px,[],1,1,3).*slocy_aa;
            aux2_3    = 2/3*mu(m)*repmat(sum(aux2_3Aux,4)./edet_aa,1,1,1,9).*slocy3b_aa;
            aux2_4Aux = reshape(py,[],1,1,3).*slocx_aa + reshape(px,[],1,1,3).*slocy_aa;
            aux2_4    = 2/3*mu(m)*repmat(sum(aux2_4Aux,4)./edet_aa,1,1,1,9).*slocx3b_aa;
            A42(:,m,1,:) = aux1_1 + aux2_1 + aux2_3;
            A42(:,m,2,:) = aux1_2 + aux2_2 + aux2_4;
            
            px = pxBackup;
            py = pyBackup;
        end
        A42      = sparse(TriInfo.ii2(:), TriInfo.jj2(:), A42(:));
        
        EdPhi = sparse(TriInfo.indicesIiPhi(:),TriInfo.indicesJjPhi(:), EdPhi(:));
        EdPhi = EdPhi + A42*u(:);
        
        % v*Theta*q -> q
        T   = zeros(nelement,1,1,9);
        for k=1:nelement
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
        T = sparse(TriInfo.ii3(:), TriInfo.jj3(:), T(:));
        
        G1 = EdPhi1*phi(:) + EdPhi - T*q;
        G2 = matrices.GradSq*phi(:);
        G3 = 0.5*matrices.Id - matrices.Mloc*phi(:);
        
        G1 = ShortenPhi(G1, TriInfo);
        G2 = ShortenPhi(G2, TriInfo);
        G3 = ShortenPhi(G3, TriInfo);
        
        G  = G1 + alpha*epsilon*G2 + alpha/epsilon*G3;
    else
        G  = [];
        G1 = [];
        G2 = [];
        G3 = [];
    end
end






function phi = ComputePhiCutoff(phi)
    
    delta = 1e-3;
    a     = 1/pi*delta^4;
    b     = pi*(1-2*delta^3)/delta^4;
    
    set1  = phi<0;
    set2  = phi>=0 & phi<=delta;
    phi(set1) = a*atan(b*phi(set1)) + delta^4;
    phi(set2) = phi(set2) - 2*delta*(phi(set2)-delta).^3 - (phi(set2)-delta).^4;
end

function cDerivative = ComputeCDerivative(phi, TriInfo)
    
    delta = 1e-3;
    a     = 1/pi*delta^4;
    b     = pi*(1-2*delta^3)/delta^4;
    
    phi   = phi(TriInfo.e2p);
    phi   = phi(:);
    set1  = phi<0;
    set2  = phi>=0 & phi <= delta;
    set3  = phi>delta;
    
    cDerivative       = zeros(size(phi));
    cDerivative(set1) = a*b ./ (1+(b.*phi(set1)).^2);
    cDerivative(set2) = 1 - 6*delta*(phi(set2)-delta).^2 - 4*(phi(set2)-delta).^3;
    cDerivative(set3) = 1;
    cDerivative       = reshape(cDerivative, [], 3);
    cDerivative       = mean(cDerivative, 2);
    cDerivative       = repmat(cDerivative, 1, 3);
end

function phi = ComputePhiCutoffEigen(phi)
    
    delta = 1e-3;
    
    set1  = phi<0;
    set2  = phi>1;
    phi(set1) = delta*atan(phi(set1)/delta);
    phi(set2) = 1 + delta*atan((phi(set2)-1)/delta);
end
