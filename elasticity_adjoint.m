function [u,p] = elasticity_adjoint(phi,Theta,TriInfo,Transformation,matrices,constants,material)
    % ELASTICITY_ADJOINT solves the elasticity and the adjoint equation for u and p
    %
    %  Input: phi            - phase-field parameter, influences elasticity tensor
    %         TriInfo        - Information of triangulation as nelement, npoint, etc
    %         Transformation - mass matrices etc of FEM
    %         matrices       - matrices of FEM needed here, see main file
    %         constants      - values of constants in a structure (gamma, M, c)
    %         material       - material properties in a structure (Lame lambda and
    %                          mu, eigenstrain eps0, thermal stress sigma0)
    %
    %  Output: u - state variable - solution of elasticity problem
    %          p - adjoint variable - solution of adjoint problem
    
    e2p            = TriInfo.e2p;
    nelement       = TriInfo.nelement;
    npoint         = TriInfo.npoint;
    idp            = TriInfo.idp;
    nphi           = TriInfo.nphi;
    sizePhi        = TriInfo.sizePhi;
    
    lambda  = material.lambda;
    mu      = material.mu;
    eps0    = material.eps0;
    sigma0  = material.sigma0;
    
    phi = ProlongPhi(phi, TriInfo);
    
    bEla = zeros(nelement,2,nphi);
    
    slocx3a_aa = matrices.slocx3a_aa;
    slocxx_aa  = matrices.slocxx_aa;
    slocy3a_aa = matrices.slocy3a_aa;
    slocyy_aa  = matrices.slocyy_aa;
    slocxy_aa  = matrices.slocxy_aa;
    slocyx_aa  = matrices.slocyx_aa;
    
    A11_2_1  = zeros(nelement,2,sizePhi,9);
    A12      = zeros(nelement,2,2,9);
    A12_11_1 = zeros(nelement,1);
    A12_11_2 = zeros(nelement,1);
    A12_22_1 = zeros(nelement,1);
    A12_22_2 = zeros(nelement,1);
    A12_12_1 = zeros(nelement,1);
    A12_12_2 = zeros(nelement,1);
    
    % C(phi)e(u):e(v) -> u
    for i=1:TriInfo.sizePhi
        phiAux = phi(:,i);
        sumPhi = sum(phiAux(e2p),2);
        A12_11_1 = A12_11_1 + (2/3*mu(i)+1/3*lambda(i))*sumPhi;
        A12_11_2 = A12_11_2 + 1/3*mu(i)*sumPhi;
        A12_22_1 = A12_22_1 + 1/3*mu(i)*sumPhi;
        A12_22_2 = A12_22_2 + (2/3*mu(i)+1/3*lambda(i))*sumPhi;
        A12_12_1 = A12_12_1 + 1/3*mu(i)*sumPhi;
        A12_12_2 = A12_12_2 + 1/3*lambda(i)*sumPhi;
    end
    A12(:,1,1,:) = repmat(A12_11_1,1,1,1,9).*slocxx_aa + repmat(A12_11_2,1,1,1,9).*slocyy_aa;
    A12(:,2,2,:) = repmat(A12_22_1,1,1,1,9).*slocxx_aa + repmat(A12_22_2,1,1,1,9).*slocyy_aa;
    A12(:,1,2,:) = repmat(A12_12_1,1,1,1,9).*slocyx_aa + repmat(A12_12_2,1,1,1,9).*slocxy_aa;
    A12(:,2,1,:) = repmat(A12_12_1,1,1,1,9).*slocxy_aa + repmat(A12_12_2,1,1,1,9).*slocyx_aa;
    % phi_i*phi1_1*tr(e(v)) -> phi, phi1
    phiAux = phi(:,1);
    phiAux = phiAux(e2p);
    for i = 1:sizePhi
        aux1x   = 1/6*eps0*(mu(i)+lambda(i))*slocx3a_aa.*reshape(phiAux(:,[1 1 1 2 2 2 3 3 3]), [], 1, 1, 9);
        aux1y   = 1/6*eps0*(mu(i)+lambda(i))*slocy3a_aa.*reshape(phiAux(:,[1 1 1 2 2 2 3 3 3]), [], 1, 1, 9);
        aux2x   = 1/6*eps0*(mu(i)+lambda(i))*slocx3a_aa.*repmat(sum(phiAux,2), 1, 1, 1, 9);
        aux2y   = 1/6*eps0*(mu(i)+lambda(i))*slocy3a_aa.*repmat(sum(phiAux,2), 1, 1, 1, 9);
        A11_2_1(:,1,i,:) = aux1x + aux2x;
        A11_2_1(:,2,i,:) = aux1y + aux2y;
    end
    
    for k=1:nelement
        clocx = Transformation{k,7};
        clocy = Transformation{k,8};
        
        phiMatrix0ele = [phi(e2p(k,1),:)',phi(e2p(k,2),:)',phi(e2p(k,3),:)'];
        
        sigmaXphi     = sigma0*phiMatrix0ele(2,:);
        bEla( k,1,:) = sigmaXphi*clocx ;
        bEla( k,2,:) = sigmaXphi*clocy ;
    end
    
    A       = sparse(TriInfo.indicesIEla(:), TriInfo.indicesJEla(:), A12(:));
    A11_2_1 = sparse(TriInfo.ii1(:), TriInfo.jj1(:), A11_2_1(:));
    bEla    = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bEla(:));
    bEla    = A11_2_1*phi(:) - bEla;
    
    
    
    
    
    
    
    
    
    
    bAdj = zeros(nelement,1,1,3);
    for k=1:nelement
        edet   = Transformation{k,1};
        slocx  = Transformation{k,9};
        slocy  = Transformation{k,10};
        mloc   = Transformation{k,11};
        
        ThetaAux = Theta(e2p(k,:));
        
        bAdj(k,1,:) = 2*slocx'/edet*(ThetaAux'*mloc*ThetaAux);
        bAdj(k,2,:) = 2*slocy'/edet*(ThetaAux'*mloc*ThetaAux);
    end
    bAdj    = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bAdj(:));
    
    
    
    
    
    
    
%     bAdj    = matrices.TrXD+matrices.TrYD;
    
    id      = ~[idp(1:npoint)==1;idp(1:npoint)==1];
    u       = zeros(2*npoint,1);
    u(id)   = A(id,id)\bEla(id);
    p       = zeros(2*npoint,1);
    p(id)   = A(id,id)\bAdj(id);
end




