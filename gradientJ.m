function [G, G1, G2, G3] = gradientJ(phi,u,p,Theta,adjoint1,TriInfo,Transformation,matrices,constants,material)
    % GRADIENTJ gives back the gradient of the objective functional J at (phi,u)
    %          using the adjoint variable p
    %
    % Input: phi            - phase-field vector
    %        u              - state variable (dependent on phi)
    %        p              - adjoint variable (dependent on u and phi)
    %        TriInfo        - Information of triangulation as nelement, npoint, etc
    %        Transformation - mass matrices etc of FEM
    %        matrices       - matrices needed here (see main file)
    %        constants      - values of constants in a structure (gamma, M, c)
    %        material       - material properties in a structure (Lame lambda and
    %                         mu, eigenstrain eps0, thermal stress sigma0)
    %
    % Output: Gradient - gradient at (phi,u)
    
    %% extract information of triangulation, material properties and constants
    e2p            = TriInfo.e2p;
    nelement       = TriInfo.nelement;
    npoint         = TriInfo.npoint;
    nphi           = TriInfo.nphi;
    sizePhi        = TriInfo.sizePhi;
    
    Mloc    = matrices.Mloc;
    Id      = matrices.Id;
    GradSq  = matrices.GradSq;
    
    alpha   = constants.alpha;
    epsilon = constants.epsilon;
    
    lambda  = material.lambda;
    mu      = material.mu;
    eps0    = material.eps0;
    sigma0  = material.sigma0;
    
    
    
    
    
    
    npoint                  = TriInfo.npoint;
    nelement                = TriInfo.nelement;
    e2p                     = TriInfo.e2p;
    id                      = ~TriInfo.idp(1:TriInfo.npoint)==1;
    sizePhi                 = TriInfo.sizePhi;
    constants.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    constants.epsilonR(2:3) = [];
    epsilonR                = constants.epsilonR;
    
    
    
    
    
    phi = ProlongPhi(phi, TriInfo);
    
    EdPhi  = zeros(nelement,sizePhi,nphi);
    
    %% Assemly matrices
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
    end
    A42      = sparse(TriInfo.ii2(:), TriInfo.jj2(:), A42(:));
    
    EdPhi = sparse(TriInfo.indicesIiPhi(:),TriInfo.indicesJjPhi(:), EdPhi(:));
    EdPhi = EdPhi + A42*u(:);
    











    % deltaPhi*Theta*q
    T = zeros(nelement,1,1,9);
    ii3      = zeros(nelement,sizePhi,1,9);
    jj3      = zeros(nelement,sizePhi,1,9);
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
        
        for i1=1:sizePhi
            for j1=1:1
                ii3( k,i1,j1,: ) = (i1-1)*npoint + ...
                    [e2p(k,1) e2p(k,2) e2p(k,3) ...
                    e2p(k,1) e2p(k,2) e2p(k,3) ...
                    e2p(k,1) e2p(k,2) e2p(k,3)];
                jj3( k,i1,j1,: ) = (j1-1)*npoint + ...
                    [e2p(k,1) e2p(k,1) e2p(k,1) ...
                    e2p(k,2) e2p(k,2) e2p(k,2) ...
                    e2p(k,3) e2p(k,3) e2p(k,3)];
            end
        end
        
        for i=1:sizePhi
            T(k,i,1,:) = 1/120*epsilonR(i)*edet*addition(:);
        end
    end
    T = sparse(ii3(:), jj3(:), T(:));
    
    gradient = -T*adjoint1;
%     gradient = reshape(gradient,[],sizePhi);
%     gradient = gradient(TriInfo.phiRowsFree,:);



% EdPhi1 = 0;
% EdPhi = 0;
% gradient = 0;


    %% And compute the gradients
    
    G1 = EdPhi1*phi(:) + EdPhi + gradient(:);
    G2 = GradSq*phi(:);
    G3 = 0.5*Id - Mloc*phi(:);
    
    G1 = ShortenPhi(G1, TriInfo);
    G2 = ShortenPhi(G2, TriInfo);
    G3 = ShortenPhi(G3, TriInfo);
    
    G = G1 + alpha*epsilon*G2 + alpha/epsilon*G3;
end