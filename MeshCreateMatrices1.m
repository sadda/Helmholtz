function [TriInfo, Transformation] = MeshCreateMatrices1(x,y,e2p,requiredX,requiredY,sizePhi,opticalCavity)
    nelement = size(e2p,1);
    npoint   = length(x);
    nphi     = 3;
    tol      = 1e-12; % For numerical errors
    
    %% Create Dirichlet boundary
    idp = zeros(size(x));
    idp(x == min(x) | x == max(x)) = 1;
    idp(y == min(y) | y == max(y)) = 1;
    
    %% Create  (required positions)
    ide       = 60*ones(size(e2p,1),1);
    ideCavity = zeros(size(e2p,1),1);
    for k=1:nelement
        % Start assigning elements (artificially prescribed air will be rewritten). AIR MUST BE LAST.
        if length(requiredX) > sizePhi && (max(x(e2p(k,:))) <= requiredX(1,sizePhi+1)+tol || min(x(e2p(k,:))) >= requiredX(2,sizePhi+1)-tol)
            ide(k) = sizePhi;
        end
        if length(requiredY) > sizePhi && (max(y(e2p(k,:))) <= requiredY(1,sizePhi+1)+tol || min(y(e2p(k,:))) >= requiredY(2,sizePhi+1)-tol)
            ide(k) = sizePhi;
        end
        for m=sizePhi:-1:1
            if min(x(e2p(k,:))) >= requiredX(1,m)-tol && max(x(e2p(k,:))) <= requiredX(2,m)+tol && min(y(e2p(k,:))) >= requiredY(1,m)-tol && max(y(e2p(k,:))) <= requiredY(2,m)+tol
                ide(k) = m;
            end
        end
        if min(x(e2p(k,:))) >= opticalCavity(1,1)-tol && max(x(e2p(k,:))) <= opticalCavity(2,1)+tol && min(y(e2p(k,:))) >= opticalCavity(1,2)-tol && max(y(e2p(k,:))) <= opticalCavity(2,2)+tol % Optical cavity
            ideCavity(k) = 1;
        end
    end
    
    %% Transformation for Triangulation
    Transformation = cell(nelement+1,11);
    
    Transformation{nelement+1,1}  = 'edet';
    Transformation{nelement+1,2}  = 'dFinv';
    Transformation{nelement+1,3}  = 'slocxx';
    Transformation{nelement+1,4}  = 'slocyy';
    Transformation{nelement+1,5}  = 'slocxy';
    Transformation{nelement+1,6}  = 'slocyx';
    Transformation{nelement+1,7}  = 'clocx';
    Transformation{nelement+1,8}  = 'clocy';
    Transformation{nelement+1,9}  = 'slocx';
    Transformation{nelement+1,10} = 'slocy';
    Transformation{nelement+1,11} = 'mloc';
    
    for k=1:nelement             % loop over elements
        [edet,dFinv]  = generatetransformation(k,e2p,x,y); % compute map
        
        [slocx,slocxx,slocy,slocyy,slocxy,mloc,clocx,clocy]=localmatrices(edet,dFinv);
        
        slocyx = slocxy';
        
        Transformation{k,1}  = edet;
        Transformation{k,2}  = dFinv;
        Transformation{k,3}  = slocxx;
        Transformation{k,4}  = slocyy;
        Transformation{k,5}  = slocxy;
        Transformation{k,6}  = slocyx;
        Transformation{k,7}  = clocx;
        Transformation{k,8}  = clocy;
        Transformation{k,9}  = slocx;
        Transformation{k,10} = slocy;
        Transformation{k,11} = mloc;
    end
    %% Generate prescribed set
    % Compute the prescribed set
    phiPrescribed = zeros(npoint,sizePhi);
    for m = sizePhi:-1:1
        for k = 1:nelement
            % Start with the less important
            if ide(k) == m
                % Reset the less important material
                phiPrescribed(e2p(k,:),:) = 0;
                % And plug in the more important
                phiPrescribed(e2p(k,:),ide(k)) = 1;
            end
        end
    end
    phiRowsFixed = sum(phiPrescribed, 2) > 0;
    phiRowsFree = ~phiRowsFixed;
    
    phiProlongationMatrix  = speye(npoint, npoint);
    phiProlongationMatrix  = phiProlongationMatrix(:,phiRowsFree);
    phiProlongationMatrix6 = kron(eye(sizePhi), phiProlongationMatrix);
    phiProlongationVector  = phiPrescribed;
    phiRowsFree6           = repmat(phiRowsFree, sizePhi, 1);
    
    npointRed              = sum(phiRowsFree);
    
    TriInfo = struct('x',x,'y',y,'requiredX',requiredX,'requiredY',requiredY,'npoint',npoint,'npointRed',npointRed,'nelement',nelement,'e2p',e2p,'idp',idp,'ide',ide,'ideCavity',ideCavity,'opticalCavity',opticalCavity,'nphi',nphi,'sizePhi',sizePhi,'phiProlongationMatrix',phiProlongationMatrix,'phiProlongationMatrix6',phiProlongationMatrix6,'phiProlongationVector',phiProlongationVector,'phiRowsFree',phiRowsFree,'phiRowsFree6',phiRowsFree6);
end

function [edet,dFinv]=generatetransformation(k,e2p,x,y)
    
    dx1=x(e2p(k,2))-x(e2p(k,1));
    dy1=y(e2p(k,2))-y(e2p(k,1));
    
    dx2=x(e2p(k,3))-x(e2p(k,1));
    dy2=y(e2p(k,3))-y(e2p(k,1));
    
    % determinant of Jacobian
    edet = dx1.*dy2 - dx2.*dy1;
    
    % inverse Jacobian
    dFinv(1,1)= dy2/edet;
    dFinv(1,2)=-dx2/edet;
    dFinv(2,1)=-dy1/edet;
    dFinv(2,2)= dx1/edet;
end

function [slocx,slocxx,slocy,slocyy,slocxy,mloc,clocx,clocy]=localmatrices(edet,dFinv)
    
    gradphip1loc(1,1) = -1;
    gradphip1loc(2,1) = -1;
    gradphip1loc(1,2) =  1;
    gradphip1loc(2,2) =  0;
    gradphip1loc(1,3) =  0;
    gradphip1loc(2,3) =  1;
    
    % dphi = ndim x nphi
    dphi = dFinv'*gradphip1loc;
    
    slocxx=1/2 * dphi(1,:)'*dphi(1,:) * edet;       % slocxx
    slocx=([1;0]'*dphi)*edet/2;                     % slocx
    slocyy=1/2 * dphi(2,:)'*dphi(2,:) * edet;       % slocyy
    slocy=([0;1]'*dphi)*edet/2;                     % slocy
    slocxy=1/2 * dphi(1,:)'*dphi(2,:) * edet;       % slocxy
    mloc=edet*[1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12;   % mloc
    clocx=[slocx;slocx;slocx]/3;                    % clocx
    clocy=[slocy;slocy;slocy]/3;                    % clocy
end
