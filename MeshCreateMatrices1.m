function [TriInfo, Transformation] = MeshCreateMatrices1(x,y,e2p,requiredX,requiredY,sizePhi)
    % Adds some information to the mesh (Dirichlet nodes, prescribed nodes), local transformations and some matrices connected to the prolongation from the free nodes
    
    nelement = size(e2p,1);
    npoint   = length(x);
    nphi     = 3;
    
    %% Dirichlet boundary
    
    idp = zeros(size(x));
    idp(x == min(x) | x == max(x)) = 1;
    idp(y == min(y) | y == max(y)) = 1;
    
    %% Generate prescribed set
    
    phiPrescribed = zeros(npoint,sizePhi);
    for i=1:npoint
        for m=1:sizePhi
            if x(i) >= requiredX(1,m) && x(i) <= requiredX(2,m) && y(i) >= requiredY(1,m) && y(i) <= requiredY(2,m)
                phiPrescribed(i,m) = 1;
            end
        end
    end
    if any(sum(phiPrescribed,2) > 1)
        error('Multiple nodes prescribed');
    end
    phiRowsFixed = sum(phiPrescribed, 2) > 0;
    phiRowsFree = ~phiRowsFixed;
    
    phiProlongationMatrix  = speye(npoint, npoint);
    phiProlongationMatrix  = phiProlongationMatrix(:,phiRowsFree);
    phiProlongationMatrix6 = kron(eye(sizePhi), phiProlongationMatrix);
    phiProlongationVector  = phiPrescribed;
    phiRowsFree6           = repmat(phiRowsFree, sizePhi, 1);
    npointRed              = sum(phiRowsFree);
    
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
    
    for k=1:nelement
        [edet,dFinv]                                        = GenerateTransformation(k,e2p,x,y);
        [slocx,slocxx,slocy,slocyy,slocxy,mloc,clocx,clocy] = LocalMatrices(edet,dFinv);
        
        Transformation{k,1}  = edet;
        Transformation{k,2}  = dFinv;
        Transformation{k,3}  = slocxx;
        Transformation{k,4}  = slocyy;
        Transformation{k,5}  = slocxy;
        Transformation{k,6}  = slocxy';
        Transformation{k,7}  = clocx;
        Transformation{k,8}  = clocy;
        Transformation{k,9}  = slocx;
        Transformation{k,10} = slocy;
        Transformation{k,11} = mloc;
    end
    
    TriInfo = struct('x',x,'y',y,'requiredX',requiredX,'requiredY',requiredY,'npoint',npoint,'npointRed',npointRed,'nelement',nelement,'e2p',e2p,'idp',idp,'nphi',nphi,'sizePhi',sizePhi,'phiProlongationMatrix',phiProlongationMatrix,'phiProlongationMatrix6',phiProlongationMatrix6,'phiProlongationVector',phiProlongationVector,'phiRowsFree',phiRowsFree,'phiRowsFree6',phiRowsFree6);
end

function [edet,dFinv] = GenerateTransformation(k,e2p,x,y)
    
    dx1   = x(e2p(k,2))-x(e2p(k,1));
    dy1   = y(e2p(k,2))-y(e2p(k,1));
    dx2   = x(e2p(k,3))-x(e2p(k,1));
    dy2   = y(e2p(k,3))-y(e2p(k,1));
    edet  = dx1.*dy2 - dx2.*dy1;           % Determinant of the Jacobian
    dFinv = 1/edet*[dy2, -dx2; -dy1, dx1]; % Inverse Jacobian
end

function [slocx,slocxx,slocy,slocyy,slocxy,mloc,clocx,clocy] = LocalMatrices(edet,dFinv)
    
    dphi   = dFinv'*[-1 1 0; -1 0 1];
    
    slocxx = 1/2 * dphi(1,:)'*dphi(1,:) * edet;       % slocxx
    slocx  = ([1;0]'*dphi)*edet/2;                    % slocx
    slocyy = 1/2 * dphi(2,:)'*dphi(2,:) * edet;       % slocyy
    slocy  = ([0;1]'*dphi)*edet/2;                    % slocy
    slocxy = 1/2 * dphi(1,:)'*dphi(2,:) * edet;       % slocxy
    mloc   = edet*[1 1/2 1/2;1/2 1 1/2;1/2 1/2 1]/12; % mloc
    clocx  = [slocx;slocx;slocx]/3;                   % clocx
    clocy  = [slocy;slocy;slocy]/3;                   % clocy
end


