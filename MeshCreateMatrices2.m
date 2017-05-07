function [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo)
    % Compute all required matrices from their local information
    
    nelement = TriInfo.nelement;
    nphi     = TriInfo.nphi;
    npoint   = TriInfo.npoint;
    e2p      = TriInfo.e2p;
    sizePhi  = TriInfo.sizePhi;
    
    %% Generate indices and empty matrices
    
    % Corresponding to sizePhi*sizePhi
    ii_Phi      = zeros(nelement,sizePhi,sizePhi,nphi^2);
    jj_Phi      = zeros(nelement,sizePhi,sizePhi,nphi^2);
    iiv_Phi     = zeros(nelement,sizePhi,nphi);
    jjv_Phi     = zeros(nelement,sizePhi,nphi);
    
    H1scal_aa   = zeros(nelement,sizePhi,sizePhi,nphi^2);
    GradSq_aa   = zeros(nelement,sizePhi,sizePhi,nphi^2);
    Mloc_aa     = zeros(nelement,sizePhi,sizePhi,nphi^2);
    Id_aa       = zeros(nelement,sizePhi,nphi);
    
    % Corresponding to sizePhi*2
    ii1         = zeros(nelement,2,sizePhi,nphi^2);
    jj1         = zeros(nelement,2,sizePhi,nphi^2);
    ii2         = zeros(nelement,sizePhi,2,nphi^2);
    jj2         = zeros(nelement,sizePhi,2,nphi^2);
    ii3         = zeros(nelement,sizePhi,1,nphi^2);
    jj3         = zeros(nelement,sizePhi,1,nphi^2);
    
    % Corresponding to 2*2 (elasticity)
    ii_ela      = zeros(nelement,2,2,nphi^2);
    jj_ela      = zeros(nelement,2,2,nphi^2);
    iiv         = zeros(nelement,2,nphi);
    jjv         = zeros(nelement,2,nphi);
    
    H1scal2D_aa = zeros(nelement,2,2,nphi^2);
    normESq_aa  = zeros(nelement,2,2,nphi^2);
    Mloc2D_aa   = zeros(nelement,2,2,nphi^2);
    Tr2D_aa     = zeros(nelement,2,2,nphi^2);
    Id2D_aa     = zeros(nelement,2,nphi);
    
    % Empty matrices for the sparse creation
    edet_aa     = zeros(nelement,1,1,1);
    slocx_aa    = zeros(nelement,1,1,nphi);
    slocx3a_aa  = zeros(nelement,1,1,nphi^2);
    slocx3b_aa  = zeros(nelement,1,1,nphi^2);
    slocxx_aa   = zeros(nelement,1,1,nphi^2);
    slocy_aa    = zeros(nelement,1,1,nphi);
    slocy3a_aa  = zeros(nelement,1,1,nphi^2);
    slocy3b_aa  = zeros(nelement,1,1,nphi^2);
    slocyy_aa   = zeros(nelement,1,1,nphi^2);
    mloc_aa     = zeros(nelement,1,1,nphi^2);
    slocxy_aa   = zeros(nelement,1,1,nphi^2);
    slocyx_aa   = zeros(nelement,1,1,nphi^2);
    slocxx1_aa  = zeros(nelement,1,1,nphi^2);
    slocxx2_aa  = zeros(nelement,1,1,nphi^2);
    slocxx3_aa  = zeros(nelement,1,1,nphi^2);
    slocyy1_aa  = zeros(nelement,1,1,nphi^2);
    slocyy2_aa  = zeros(nelement,1,1,nphi^2);
    slocyy3_aa  = zeros(nelement,1,1,nphi^2);
    slocxy1_aa  = zeros(nelement,1,1,nphi^2);
    slocxy2_aa  = zeros(nelement,1,1,nphi^2);
    slocxy3_aa  = zeros(nelement,1,1,nphi^2);
    slocyx1_aa  = zeros(nelement,1,1,nphi^2);
    slocyx2_aa  = zeros(nelement,1,1,nphi^2);
    slocyx3_aa  = zeros(nelement,1,1,nphi^2);
    
    for k=1:nelement
        edet   = Transformation{k,1};
        slocxx = Transformation{k,3};
        slocyy = Transformation{k,4};
        slocxy = Transformation{k,5};
        slocyx = Transformation{k,6};
        clocx  = Transformation{k,7};
        clocy  = Transformation{k,8};
        slocx  = Transformation{k,9};
        slocy  = Transformation{k,10};
        mloc   = Transformation{k,11};
        
        %% Fill indices
        
        e2pRow1 = repmat(e2p(k,:), 1, 3);
        e2pRow2 = kron(e2p(k,:), ones(1, 3));
        for m=1:sizePhi
            for j1=1:sizePhi
                ii_Phi(k,m,j1,:) = (m-1)*npoint + e2pRow1;
                jj_Phi(k,m,j1,:) = (j1-1)*npoint + e2pRow2;
            end
            iiv_Phi(k,m,:) = (m-1)*npoint + e2p(k,1:3);
            jjv_Phi(k,m,:) = ones(1,3);
        end
        for i1=1:2
            for j1=1:sizePhi
                ii1(k,i1,j1,:) = (i1-1)*npoint + e2pRow1;
                jj1(k,i1,j1,:) = (j1-1)*npoint + e2pRow2;
            end
        end
        for i1=1:sizePhi
            for j1=1:2
                ii2(k,i1,j1,:) = (i1-1)*npoint + e2pRow1;
                jj2(k,i1,j1,:) = (j1-1)*npoint + e2pRow2;
            end
        end
        for i1=1:sizePhi
            for j1=1:1
                ii3(k,i1,j1,:) = (i1-1)*npoint + e2pRow1;
                jj3(k,i1,j1,:) = (j1-1)*npoint + e2pRow2;
            end
        end
        for i1=1:2
            for j1=1:2
                ii_ela(k,i1,j1,:) = (i1-1)*npoint + e2pRow1;
                jj_ela(k,i1,j1,:) = (j1-1)*npoint + e2pRow2;
            end
        end
        iiv(k,1,:) =          e2p(k,1:3);
        iiv(k,2,:) = npoint + e2p(k,1:3);
        jjv(k,1,:) = ones(1,3);
        jjv(k,2,:) = ones(1,3);
        
        %% Fill matrices
        
        % If there are two indices (like m and j1), it corresponds to a matrix. If there is only one index (like m), it corresponds to a vector.
        % Phi is represented in matrix (npoint,4) and when a vector is created, the first npoint indices will correspond to germanium. The same holds true for u (first npoint displacements along x axis).
        % For ii_Phi, the second and third index corresponds to all (4*4) possible combination of two different materials. The last index is the number of nodes corresponding to the given element.
        % For Mloc2D_aa and Tr2D_aa we want to consider only (.,1,1,.) and (.,2,2,.) because there is no combination of u_x and v_y.
        
        for m=1:sizePhi
            H1scal_aa(k,m,m,:) = mloc(:)+slocxx(:)+slocyy(:);
            GradSq_aa(k,m,m,:) = slocxx(:)+slocyy(:);
            Mloc_aa(k,m,m,:)   = mloc(:);
            Id_aa(k,m,:)       = 1/6*edet*ones(3,1);
        end
        H1scal2D_aa(k,1,1,:) = mloc(:)+slocxx(:)+slocyy(:);
        H1scal2D_aa(k,2,2,:) = mloc(:)+slocxx(:)+slocyy(:);
        normESq_aa(k,1,1,:)  = slocxx(:)+0.5*slocyy(:);
        normESq_aa(k,1,2,:)  = 0.5*slocyx(:);
        normESq_aa(k,2,2,:)  = slocyy(:)+0.5*slocxx(:);
        normESq_aa(k,2,1,:)  = 0.5*slocxy(:);
        Mloc2D_aa(k,1,1,:)   = mloc(:);
        Mloc2D_aa(k,2,2,:)   = mloc(:);
        Tr2D_aa(k,1,1,:)     = clocx(:);
        Tr2D_aa(k,2,2,:)     = clocy(:);
        Id2D_aa(k,1,:)       = 1/6*edet*ones(3,1);
        Id2D_aa(k,2,:)       = 1/6*edet*ones(3,1);
        
        edet_aa(k,1,1,1)     = edet;
        mloc_aa(k,1,1,:)     = mloc(:);
        slocx_aa(k,1,1,:)    = slocx(:);
        slocx3a_aa(k,1,1,:)  = repmat(slocx(:), 3, 1);
        slocx3b_aa(k,1,1,:)  = reshape(repmat(slocx, 3, 1), 9 , 1);
        slocxx_aa(k,1,1,:)   = slocxx(:);
        slocxx1_aa(k,1,1,:)  = repmat(slocxx(:,1), 3, 1);
        slocxx2_aa(k,1,1,:)  = repmat(slocxx(:,2), 3, 1);
        slocxx3_aa(k,1,1,:)  = repmat(slocxx(:,3), 3, 1);
        slocy_aa(k,1,1,:)    = slocy(:);
        slocy3a_aa(k,1,1,:)  = repmat(slocy(:), 3, 1);
        slocy3b_aa(k,1,1,:)  = reshape(repmat(slocy, 3, 1), 9, 1);
        slocyy_aa(k,1,1,:)   = slocyy(:);
        slocyy1_aa(k,1,1,:)  = repmat(slocyy(:,1), 3, 1);
        slocyy2_aa(k,1,1,:)  = repmat(slocyy(:,2), 3, 1);
        slocyy3_aa(k,1,1,:)  = repmat(slocyy(:,3), 3, 1);
        slocxy_aa(k,1,1,:)   = slocxy(:);
        slocxy1_aa(k,1,1,:)  = repmat(slocxy(:,1), 3, 1);
        slocxy2_aa(k,1,1,:)  = repmat(slocxy(:,2), 3, 1);
        slocxy3_aa(k,1,1,:)  = repmat(slocxy(:,3), 3, 1);
        slocyx_aa(k,1,1,:)   = slocyx(:);
        slocyx1_aa(k,1,1,:)  = repmat(slocyx(:,1), 3, 1);
        slocyx2_aa(k,1,1,:)  = repmat(slocyx(:,2), 3, 1);
        slocyx3_aa(k,1,1,:)  = repmat(slocyx(:,3), 3, 1);
    end
    
    %% Create matrices from sparse data
    
    H1scal   = sparse(ii_Phi(:),jj_Phi(:),H1scal_aa(:));
    GradSq   = sparse(ii_Phi(:),jj_Phi(:),GradSq_aa(:));
    Mloc     = sparse(ii_Phi(:),jj_Phi(:),Mloc_aa(:));
    Id       = sparse(iiv_Phi(:),jjv_Phi(:),Id_aa(:));
    H1scal2D = sparse(ii_ela(:),jj_ela(:),H1scal2D_aa(:));
    normESq  = sparse(ii_ela(:),jj_ela(:),normESq_aa(:));
    Mloc2D   = sparse(ii_ela(:),jj_ela(:),Mloc2D_aa(:));
    Tr2D     = sparse(ii_ela(:),jj_ela(:),Tr2D_aa(:));
    Id2D     = sparse(iiv(:),jjv(:),Id2D_aa(:));
    
    %% Save everything
    
    matrices = struct('Mloc',Mloc,'Id',Id,'H1scal',H1scal,'H1scal2D',H1scal2D,...
        'Mloc2D',Mloc2D,'Tr2D',Tr2D,'normESq',normESq,'Id2D',Id2D,'GradSq',GradSq,...
        'edet_aa',edet_aa,'slocx_aa',slocx_aa,'slocx3a_aa',slocx3a_aa,'slocx3b_aa',slocx3b_aa,'slocxx_aa',slocxx_aa,'slocy_aa',slocy_aa,'slocy3a_aa',slocy3a_aa,'slocy3b_aa',slocy3b_aa,'slocyy_aa',slocyy_aa,'mloc_aa',mloc_aa,'slocxy_aa',slocxy_aa,'slocyx_aa',slocyx_aa,...
        'slocxx1_aa',slocxx1_aa,'slocxx2_aa',slocxx2_aa,'slocxx3_aa',slocxx3_aa,'slocyy1_aa',slocyy1_aa,'slocyy2_aa',slocyy2_aa,'slocyy3_aa',slocyy3_aa,'slocxy1_aa',slocxy1_aa,'slocxy2_aa',slocxy2_aa,'slocxy3_aa',slocxy3_aa,'slocyx1_aa',slocyx1_aa,'slocyx2_aa',slocyx2_aa,'slocyx3_aa',slocyx3_aa);
    
    TriInfo.indicesIiPhi = iiv_Phi;
    TriInfo.indicesJjPhi = jjv_Phi;
    TriInfo.indicesIPhi  = ii_Phi;
    TriInfo.indicesJPhi  = jj_Phi;
    TriInfo.ii1          = ii1;
    TriInfo.ii2          = ii2;
    TriInfo.ii3          = ii3;
    TriInfo.jj1          = jj1;
    TriInfo.jj2          = jj2;
    TriInfo.jj3          = jj3;
    TriInfo.indicesIEla  = ii_ela;
    TriInfo.indicesJEla  = jj_ela;
    TriInfo.indicesIElav = iiv;
    TriInfo.indicesJElav = jjv;
end

