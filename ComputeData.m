function [J, G, J1, J2, J3, G1, G2, G3, u, Theta, dataEigen] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, dataEigen)
    
    if isempty(dataEigen)
        dataEigen = containers.Map('KeyType','double','ValueType','any');
        dataEigen(-1) = 0;
    end
    
    options  = SetField(options, 'computeG', 1);
    options  = SetField(options, 'symmetrize', 0);
    options  = SetField(options, 'jointObjectiveThetaLp', 2);
    options  = SetField(options, 'jointObjectivePhiLp', 0);
    options  = SetField(options, 'normalizationLp', 2);
    options  = SetField(options, 'cutoffs', 1);
    options  = SetField(options, 'alphaGe', 0);
    options  = SetField(options, 'geLevel2', 0);
    options  = SetField(options, 'regThetaL1', 0);
    options  = SetField(options, 'regPhiL1', 0);
    options  = SetField(options, 'regPhiL2', 0);
    
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    npoint   = TriInfo.npoint;
    idp      = TriInfo.idp;
    sizePhi  = TriInfo.sizePhi;
    epsilonR = material.epsilonR;
    
    id             = ~[idp(1:npoint)==1;idp(1:npoint)==1];
    id1D           = ~[idp(1:npoint)==1];
    phi            = ProlongPhi(phi, TriInfo);
    
    %% Compute u
    if options.method == 0 || options.method == 1
        lambda   = material.lambda;
        mu       = material.mu;
        eps0     = material.eps0;
        sigma0   = material.sigma0;
        
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
            phiAux = ComputePhiCutoff(phiAux, options);
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
        if options.originalEps0
            phiAux = phi(:,1);
        else
            phiAux = phi(:,3);
        end
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
            bEla(k,1,:)   = sigmaXphi*clocx ;
            bEla(k,2,:)   = sigmaXphi*clocy ;
        end
        
        AEla       = sparse(TriInfo.indicesIEla(:), TriInfo.indicesJEla(:), AEla12(:));
        AEla11_2_1 = sparse(TriInfo.ii1(:), TriInfo.jj1(:), AEla11_2_1(:));
        bEla       = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bEla(:));
        bEla       = AEla11_2_1*phi(:) - bEla;
        
        if options.originalEps0
            u       = zeros(2*npoint,1);
            u(id)   = AEla(id,id)\bEla(id);
        else
            uFixed  = material.eps0*[TriInfo.x; TriInfo.y];
            bEla    = bEla - AEla(:,~id)*uFixed(~id);
            
            u       = zeros(2*npoint,1);
            u(id)   = AEla(id,id)\bEla(id);
            u(~id)  = uFixed(~id);
        end
    else
        u = [];
    end
    
    %% Compute Theta and eigen
    
    if options.method == 0 || options.method == 2
        phiSum = ComputePhiCutoffEigen(phi, options);
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
        TFull = sparse(ii1D(:), jj2D(:), T(:));
        T     = TFull(id1D,id1D);
        S     = matrices.GradSq(id1D,id1D);
        M1    = matrices.Mloc(1:npoint,1:npoint);
        M     = matrices.Mloc(id1D,id1D);
        
        constM = 1+1e-3;
        if isKey(dataEigen, dataEigen(-1))
            data        = dataEigen(dataEigen(-1));
            lambdaEigen = data{1};
            phiEigen    = data{2};
            constL      = max(epsilonR);
            constCP     = 2*pi^2/((max(TriInfo.x)-min(TriInfo.x))*(max(TriInfo.y)-min(TriInfo.y)));
            constEigen  = constL*(constCP^2+1)*(2*constM*constCP^2+1)/(constCP^2);
            phiDiff     = phi - phiEigen;
            shift2      = -lambdaEigen + constEigen*sqrt(phiDiff(:)'*matrices.Mloc*phiDiff(:));
            shift2      = shift2 + 1e-5;
        else
            shift2      = Inf;
        end
        
        shift1  = constM*sum(epsilonR);
        shift   = min(shift1, shift2);
        
        [x, eigen] = eigs(S-T+shift*M, M, 1, 'sm');
        
        Theta       = zeros(npoint,1);
        Theta(id1D) = x;
        eigen       = eigen - shift;
        
        normCon1    = options.normalizationLp;
        normCon2    = options.jointObjectiveThetaLp;
        normCon3    = options.jointObjectivePhiLp;
        
        switch options.normalizationLp
            case 1
                ThetaNorm = ones(npoint,1)'*M1*Theta;
            case 2
                ThetaNorm = sqrt(Theta'*M1*Theta);
            otherwise
                ThetaNorm = Theta'.^(normCon1/2)*M1*Theta.^(normCon1/2);
                ThetaNorm = ThetaNorm^(1/normCon1);
        end
        Theta = Theta / ThetaNorm;
        if mean(Theta) < 0
            Theta = -Theta;
        end
        
        if options.symmetrize
            Theta  = SymmetryCompute(Theta, TriInfo, 0, 0, 1e-6);
        end
        dataEigen(dataEigen(-1)) = {eigen, phi};
    end
    
    %% Compute objective
    
    if options.method == 2
        J1  = 0.5*(Theta-TriInfo.phiGe)'*M1*(Theta-TriInfo.phiGe) + options.regThetaL1*ones(npoint,1)'*M1*Theta;
        J2  = ones(npoint,1)'*M1*phi(:,1);
        J3  = 0.5*phi(:,1)'*M1*phi(:,1);
        
        J   = J1 + options.regPhiL1*J2 + options.regPhiL2*J3;
    elseif options.method == 1
        
        
        % TODO dopsat
        
        blabla
        
        
    elseif options.method == 0
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
        
        if normCon3 > 0
            J1       = -Theta'.^normCon2*AObj*phi(:,1).^normCon3;
        else
            switch normCon2
                case 1
                    J1       = -ones(npoint,1)'*AObj*Theta;
                case 2
                    J1       = -Theta'*AObj*Theta;
                otherwise
                    J1       = -Theta'.^(normCon2/2)*AObj*Theta.^(normCon2/2);
            end
        end
        J2   = 0.5*(matrices.GradSq*phi(:))'*phi(:);
        J3   = 0.5*(matrices.Id'*phi(:) - (matrices.Mloc*phi(:))'*phi(:));
        
        J    = J1 + constants.alpha*constants.epsilon*J2 + constants.alpha/constants.epsilon*J3;
        
        if options.alphaGe > 0
            J = J + options.alphaGe * 0.5*(matrices.Id(1:npoint)'*phi(:,1) - phi(:,1)'*M1*phi(:,1));
        end
    end
    if options.geLevel2 > 0
        J1Add = options.gammaPen*0.5*(ones(1,npoint)*M1*phi(:,1) - options.geLevel2)^2;
        J1    = J1 + J1Add;
        J     = J + J1Add;
    end
    
    if nargout > 1 && options.computeG
        %% Compute q
        
        if options.method == 0 || options.method == 2
            % v*Theta*q -> q
            T2   = zeros(nelement,1,1,9);
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
                    T2(k,i,1,:) = 1/120*epsilonR(i)*edet*addition(:);
                end
            end
            T2  = sparse(TriInfo.ii3(:), TriInfo.jj3(:), T2(:));
            
            B11 = S-T-eigen*M;
            B12 = sparse(sum(id1D),sum(~id1D));
            B13 = -M*Theta(id1D);
            B21 = B12';
            B22 = speye(sum(~id1D));
            B23 = sparse(sum(~id1D),1);
            switch options.normalizationLp
                case 1
                    B31 = ones(1,npoint)*M1(:,id1D);
                    B32 = ones(1,npoint)*M1(:,~id1D);
                case 2
                    B31 = B13';
                    B32 = B23';
                otherwise
                    BP  = Theta.^(normCon1/2)'*M1*sparse(1:npoint,1:npoint,Theta.^(normCon1/2-1));
                    B31 = BP(:,id1D);
                    B32 = BP(:,~id1D);
            end
            B33 = 0;
            B   = [B11 B12 B13; B21 B22 B23; B31 B32 B33];
            
            if options.method == 2
                bAdjQ1   = -M1(id1D,:)*(Theta-TriInfo.phiGe) - options.regThetaL1*M1(id1D,:)*ones(npoint,1);
                bAdjQ2   = -M1(~id1D,:)*(Theta-TriInfo.phiGe) - options.regThetaL1*M1(~id1D,:)*ones(npoint,1);
                
                qAux     = B'\[bAdjQ1; bAdjQ2; 0];
                qAux     = qAux(1:end-1);
                q        = zeros(npoint,1);
                q(id1D)  = qAux(1:sum(id1D));
                q(~id1D) = qAux(sum(id1D)+1:end);
            elseif options.method == 0
                if normCon3 > 0
                    bAdjQ  = -normCon2*(Theta.^(normCon2-1)).*(AObj*phi(:,1).^normCon3);
                    bAdjQ1 = bAdjQ(id1D);
                    bAdjQ2 = bAdjQ(~id1D);
                else
                    switch normCon2
                        case 1
                            bAdjQ1 = -AObj(id1D,:)*ones(npoint,1);
                            bAdjQ2 = -AObj(~id1D,:)*ones(npoint,1);
                        case 2
                            bAdjQ1 = -2*AObj(id1D,:)*Theta;
                            bAdjQ2 = -2*AObj(~id1D,:)*Theta;
                        otherwise
                            bAdjQP = -normCon2*sparse(1:npoint,1:npoint,Theta.^(normCon2/2-1))*AObj*Theta.^(normCon2/2);
                            bAdjQ1 = bAdjQP(id1D);
                            bAdjQ2 = bAdjQP(~id1D);
                    end
                end
                bAdjQ1 = -bAdjQ1;
                bAdjQ2 = -bAdjQ2;
                
                qAux     = B'\[bAdjQ1; bAdjQ2; 0];
                qAux     = qAux(1:end-1);
                q        = zeros(npoint,1);
                q(id1D)  = qAux(1:sum(id1D));
                q(~id1D) = qAux(sum(id1D)+1:end);
            end
        end
        
        %% Compute p
        
        if options.method == 1 || options.method == 0
            if options.method == 1
                
                blabla
                
            elseif options.method == 0
                %% Compute p
                
                % tr(e(v))theta^2 -> v (vector)
                bAdjP = zeros(nelement,1,3);
                for k=1:nelement
                    edet   = Transformation{k,1};
                    slocx  = Transformation{k,9};
                    slocy  = Transformation{k,10};
                    mloc   = Transformation{k,11};
                    
                    if normCon3 > 0
                        phiAux       = phi(:,1);
                        phiAux       = phiAux(e2p(k,:));
                        ThetaAux     = Theta(e2p(k,:));
                        bAdjP(k,1,:) = 2*slocx'/edet*((phiAux.^normCon3)'*mloc*ThetaAux.^normCon2);
                        bAdjP(k,2,:) = 2*slocy'/edet*((phiAux.^normCon3)'*mloc*ThetaAux.^normCon2);
                    else
                        switch normCon2
                            case 1
                                ThetaAux     = Theta(e2p(k,:));
                                bAdjP(k,1,:) = 2*slocx'/edet*(ones(size(ThetaAux))'*mloc*ThetaAux);
                                bAdjP(k,2,:) = 2*slocy'/edet*(ones(size(ThetaAux))'*mloc*ThetaAux);
                            case 2
                                ThetaAux     = Theta(e2p(k,:));
                                bAdjP(k,1,:) = 2*slocx'/edet*(ThetaAux'*mloc*ThetaAux);
                                bAdjP(k,2,:) = 2*slocy'/edet*(ThetaAux'*mloc*ThetaAux);
                            otherwise
                                ThetaAux     = Theta(e2p(k,:)).^(normCon2/2);
                                bAdjP(k,1,:) = 2*slocx'/edet*(ThetaAux'*mloc*ThetaAux);
                                bAdjP(k,2,:) = 2*slocy'/edet*(ThetaAux'*mloc*ThetaAux);
                        end
                    end
                end
                bAdjP   = sparse(TriInfo.indicesIElav(:), TriInfo.indicesJElav(:), bAdjP(:));
                
                p       = zeros(2*npoint,1);
                p(id)   = AEla(id,id)\bAdjP(id);
            end
        end
        
        %% Compute gradient
        
        if options.method == 2
            G1 = -T2(:,id1D)*q(id1D);
            G2 = zeros(size(G1));
            G3 = zeros(size(G2));
            G2(1:npoint) = M1*ones(npoint,1);
            G3(1:npoint) = M1*phi(:,1);
            
            G1 = ShortenPhi(G1, TriInfo);
            G2 = ShortenPhi(G2, TriInfo);
            G3 = ShortenPhi(G3, TriInfo);
            
            G  = G1 + options.regPhiL1*G2 + options.regPhiL2*G3;
        elseif options.method == 1
            
            blabla
            
        elseif options.method == 0
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
            if options.originalEps0
                aux1 = 8*(lambda(1)+mu(1))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
                A41_4(:,1,1,:) = repmat(aux1,1,1,1,9).*mloc_aa;
                for m=[2 3 4]
                    aux2 = 4*(lambda(m)+mu(m))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
                    A41_4(:,1,m,:) = repmat(aux2,1,1,1,9).*mloc_aa;
                    A41_4(:,m,1,:) = repmat(aux2,1,1,1,9).*mloc_aa;
                end
            else
                aux1 = 8*(lambda(3)+mu(3))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
                A41_4(:,3,3,:) = repmat(aux1,1,1,1,9).*mloc_aa;
                for m=[1 2 4]
                    aux2 = 4*(lambda(m)+mu(m))*eps0 * sum(reshape(px,[],1,1,3).*slocx_aa+reshape(py,[],1,1,3).*slocy_aa,4) ./ edet_aa;
                    A41_4(:,3,m,:) = repmat(aux2,1,1,1,9).*mloc_aa;
                    A41_4(:,m,3,:) = repmat(aux2,1,1,1,9).*mloc_aa;
                end
            end
            EdPhi1    = sparse(TriInfo.indicesIPhi(:), TriInfo.indicesJPhi(:), -A41_4(:));
            % C(w)e(u):e(p) -> u
            for m=1:sizePhi
                pxBackup    = px;
                pyBackup    = py;
                cDerivative = ComputeCDerivative(phi(:,m), TriInfo, options);
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
            
            G1 = EdPhi1*phi(:) + EdPhi - T2(:,id1D)*q(id1D);
            G2 = matrices.GradSq*phi(:);
            G3 = 0.5*matrices.Id - matrices.Mloc*phi(:);
            if normCon3 > 0
                G1(1:npoint) = G1(1:npoint) - normCon3*(phi(:,1).^(normCon3-1)).*(AObj*Theta.^normCon2);
            end
            
            G1 = ShortenPhi(G1, TriInfo);
            G2 = ShortenPhi(G2, TriInfo);
            G3 = ShortenPhi(G3, TriInfo);
            
            G  = G1 + constants.alpha*constants.epsilon*G2 + constants.alpha/constants.epsilon*G3;
            
            if options.alphaGe > 0
                GAux = options.alphaGe*(0.5*matrices.Id(1:npoint) - M1*phi(:,1));
                GAux = TriInfo.phiProlongationMatrix'*GAux;
                G(1:TriInfo.npointRed) = G(1:TriInfo.npointRed) + GAux;
            end
        end
        if options.geLevel2 > 0
            G1Add           = zeros(sizePhi*npoint,1);
            G1Add(1:npoint) = options.gammaPen*(ones(1,npoint)*M1*phi(:,1) - options.geLevel2)*M1*ones(npoint,1);
            G1Add           = ShortenPhi(G1Add, TriInfo);
            G1(1:npoint)    = G1(1:npoint) + G1Add(1:npoint);
            G(1:npoint)     = G(1:npoint) + G1Add(1:npoint);
        end
    else
        G  = [];
        G1 = [];
        G2 = [];
        G3 = [];
    end
end







function phi = ComputePhiCutoff(phi, options)
    
    if options.cutoffs
        delta = 1e-3;
        
        a     = 1/pi*delta^4;
        b     = pi*(1-2*delta^3)/delta^4;
        
        set1  = phi<0;
        set2  = phi>=0 & phi<=delta;
        phi(set1) = a*atan(b*phi(set1)) + delta^4;
        phi(set2) = phi(set2) - 2*delta*(phi(set2)-delta).^3 - (phi(set2)-delta).^4;
    end
end

function cDerivative = ComputeCDerivative(phi, TriInfo, options)
    
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
    
    if ~options.cutoffs
        cDerivative = ones(size(cDerivative));
    end
end

function phi = ComputePhiCutoffEigen(phi, options)
    
    if options.cutoffs
        delta = 1e-3;
        
        set1  = phi<0;
        set2  = phi>1;
        phi(set1) = delta*atan(phi(set1)/delta);
        phi(set2) = 1 + delta*atan((phi(set2)-1)/delta);
    end
end
