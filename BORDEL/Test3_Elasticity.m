% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER
% NOT MAINTAINED. USE AT YOUR OWN DANGER



clear all;

%% Data preparation
CreateMatrices;

testNumber = 3;
testPhiSimple = 0;
sigma0 = 2;
eps0 = 3;
mu = 2*ones(1, 6);
lambda = 4*ones(1, 6);

material = struct('lambda', lambda,'mu',mu,'eps0',eps0,'sigma0',sigma0);

x              = TriInfo.x;
y              = TriInfo.y;
e2p            = TriInfo.e2p;
nelement       = TriInfo.nelement;
npoint         = TriInfo.npoint;
idp            = TriInfo.idp;
nphi           = TriInfo.nphi;

%% Symbolic computations

phi = zeros(npoint, 6);
if testPhiSimple == 1
    phi(:,1) = (x+2)/4;
    phi(:,4) = (2-x)/4;
else
    phi(:,1) = (x+2)/10 + y/10;
    phi(:,4) = (2-x)/20 + 2/10 - y/20;
    phi(:,2) = 1 - sum(phi,2);
end
s0 = sigma0;
e0 = eps0;
l = lambda(1);
m = mu(1);
if testNumber == 1
    % Test 1 for Dirichlet boundary only
    uExpected1 = sin(pi*x).*sin(pi/3*y);
    uExpected2 = sin(pi*x).*sin(pi/3*y);
    uExpected = [uExpected1; uExpected2];
    if testPhiSimple == 1
        fMod1 = sigma0/4 + 2*mu(1)*(eps0/4 + pi^2*sin(pi*x).*sin((pi*y)/3)) + lambda(1)*(eps0/2 - (pi^2*cos(pi*x).*cos((pi*y)/3))/3 + pi^2*sin(pi*x).*sin((pi*y)/3)) - 2*mu(1)*((pi^2*cos(pi*x).*cos((pi*y)/3))/6 - (pi^2*sin(pi*x).*sin((pi*y)/3))/18);
        fMod2 = (2*pi^2*mu(1)*sin(pi*x).*sin((pi*y)/3))/9 - 2*mu(1)*((pi^2*cos(pi*x).*cos((pi*y)/3))/6 - (pi^2*sin(pi*x).*sin((pi*y)/3))/2) - lambda(1)*((pi^2*cos(pi*x).*cos((pi*y)/3))/3 - (pi^2*sin(pi*x).*sin((pi*y)/3))/9);
    else
        fMod1 = s0/20 + 2.*m.*(e0/10 + pi^2.*sin(pi.*x).*sin((pi.*y)/3)) + l.*(e0/5 - (pi^2.*cos(pi.*x).*cos((pi.*y)/3))/3 + pi^2.*sin(pi.*x).*sin((pi.*y)/3)) - 2.*m.*((pi^2.*cos(pi.*x).*cos((pi.*y)/3))/6 - (pi^2.*sin(pi.*x).*sin((pi.*y)/3))/18);
        fMod2 = s0/20 + 2.*m.*(e0/10 + (pi^2.*sin(pi.*x).*sin((pi.*y)/3))/9) + l.*(e0/5 - (pi^2.*cos(pi.*x).*cos((pi.*y)/3))/3 + (pi^2.*sin(pi.*x).*sin((pi.*y)/3))/9) - 2.*m.*((pi^2.*cos(pi.*x).*cos((pi.*y)/3))/6 - (pi^2.*sin(pi.*x).*sin((pi.*y)/3))/2);
    end
    fMod = [fMod1; fMod2];
elseif testNumber == 2
    % Test 2 for Dirichlet boundary only
    error('Not implemeted');
    ux = (xx.^2-4)*((yy-3/2).^2-9/4);
    uy = (xx.^2-4)*((yy-3/2).^2-9/4);
else
    % Tests for Dirichlet and Neumann boundary
    uExpected1 = (sin(pi*x)).^2.*(sin(pi/3*y)).^2;
    uExpected2 = (sin(pi*x)).^2.*(sin(pi/3*y)).^2;
    uExpected = [uExpected1; uExpected2];
    if testPhiSimple == 1
        fMod1 = s0/4 + l*(e0/2 + 2*pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2 - 2*pi^2*cos(pi*x).^2.*sin((pi*y)/3).^2 - (4*pi^2*cos(pi*x).*cos((pi*y)/3).*sin(pi*x).*sin((pi*y)/3))/3) - 2*m*((pi^2*cos((pi*y)/3).^2.*sin(pi*x).^2)/9 - (pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2)/9 + (2*pi^2*cos(pi*x).*cos((pi*y)/3).*sin(pi*x).*sin((pi*y)/3))/3) + 2*m*(e0/4 + 2*pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2 - 2*pi^2*cos(pi*x).^2.*sin((pi*y)/3).^2);
        fMod2 = 2*m*((2*pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2)/9 - (2*pi^2*cos((pi*y)/3).^2.*sin(pi*x).^2)/9) - l*((2*pi^2*cos((pi*y)/3).^2.*sin(pi*x).^2)/9 - (2*pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2)/9 + (4*pi^2*cos(pi*x).*cos((pi*y)/3).*sin(pi*x).*sin((pi*y)/3))/3) - 2*m*(pi^2*cos(pi*x).^2.*sin((pi*y)/3).^2 - pi^2*sin(pi*x).^2.*sin((pi*y)/3).^2 + (2*pi^2*cos(pi*x).*cos((pi*y)/3).*sin(pi*x).*sin((pi*y)/3))/3);
    else
        fMod1 = s0/20 + l.*(e0/5 + 2.*pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2 - 2.*pi.^2.*cos(pi.*x).^2.*sin((pi.*y)/3).^2 - (4.*pi.^2.*cos(pi.*x).*cos((pi.*y)/3).*sin(pi.*x).*sin((pi.*y)/3))/3) - 2.*m.*((pi.^2.*cos((pi.*y)/3).^2.*sin(pi.*x).^2)/9 - (pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2)/9 + (2.*pi.^2.*cos(pi.*x).*cos((pi.*y)/3).*sin(pi.*x).*sin((pi.*y)/3))/3) + 2.*m.*(e0/10 + 2.*pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2 - 2.*pi.^2.*cos(pi.*x).^2.*sin((pi.*y)/3).^2);
        fMod2 = s0/20 + l.*(e0/5 + (2.*pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2)/9 - (2.*pi.^2.*cos((pi.*y)/3).^2.*sin(pi.*x).^2)/9 - (4.*pi.^2.*cos(pi.*x).*cos((pi.*y)/3).*sin(pi.*x).*sin((pi.*y)/3))/3) - 2.*m.*(pi.^2.*cos(pi.*x).^2.*sin((pi.*y)/3).^2 - pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2 + (2.*pi.^2.*cos(pi.*x).*cos((pi.*y)/3).*sin(pi.*x).*sin((pi.*y)/3))/3) + 2.*m.*(e0/10 + (2.*pi.^2.*sin(pi.*x).^2.*sin((pi.*y)/3).^2)/9 - (2.*pi.^2.*cos((pi.*y)/3).^2.*sin(pi.*x).^2)/9);
    end
    fMod = [fMod1; fMod2];
end



% syms xx yy l m e0 s0;
% ux = sin(pi*xx).*sin(pi/3*yy);
% uy = sin(pi*xx).*sin(pi/3*yy);
% % ux = (sin(pi*xx)).^2.*(sin(pi/3*yy)).^2;
% % uy = (sin(pi*xx)).^2.*(sin(pi/3*yy)).^2;
% phi1Function = (xx+2)/10 + yy/10;
% phi4Function = (2-xx)/20 + 2/10 - yy/20;
% phi = zeros(npoint, 6);
% phi(:,1) = double(subs(phi1Function, {'xx', 'yy'}, {x, y}));
% phi(:,4) = double(subs(phi4Function, {'xx', 'yy'}, {x, y}));
% eu = [diff(ux,'xx') diff(ux,'yy'); diff(uy,'xx') diff(uy,'yy')];
% eu = 1/2*(eu+transpose(eu));
% euShift = eu - e0*[phi1Function 0; 0 phi1Function];
% midPart = l*trace(euShift)*eye(2) + 2*m*euShift + s0*[phi4Function 0; 0 phi4Function];
% div = [diff(midPart(1,1),'xx') + diff(midPart(2,1),'yy'); diff(midPart(1,2),'xx') + diff(midPart(2,2),'yy')];
% f = subs(-div, {'l' 'm' 'e0', 's0'}, {material.lambda(1) material.mu(1) material.eps0 material.sigma0});
% fModCheck = [double(subs(f(1), {'xx' 'yy'}, {x y})); double(subs(f(2), {'xx' 'yy'}, {x y}))];



%% Solve Elasticity

ii = TriInfo.indicesIEla; % indices local-to-global
jj = TriInfo.indicesJEla; % indices local-to-global

aa = zeros(nelement,2,2,nphi^2); % entry of bilinear form

iiv = zeros(nelement,2,nphi); % sparse i-index
jjv = zeros(nelement,2,nphi); % sparse j-index
bElav1 = zeros(nelement,2,nphi); % entry of linear form
bElavF = zeros(size(bElav1));
for k=1:nelement             % loop over elements
    edet = Transformation{k,1};
    slocxx = Transformation{k,3};
    slocyy = Transformation{k,4};
    slocxy = Transformation{k,5};
    slocyx = Transformation{k,6};
    clocx = Transformation{k,7};
    clocy = Transformation{k,8};
    slocx = Transformation{k,9};
    slocy = Transformation{k,10};
    
    iiv( k,1,:) =          e2p(k,1:3);
    iiv( k,2,:) = npoint + e2p(k,1:3);
    jjv( k,1,:) = ones(1,3);
    jjv( k,2,:) = ones(1,3);
    
    % Lame coefficients on element
    phiMatrix0ele = [phi(e2p(k,1),:)',phi(e2p(k,2),:)',phi(e2p(k,3),:)']; % coefficients of phi on the element
    
    muXphi = mu*phiMatrix0ele;
    laXphi = lambda*phiMatrix0ele;
    
    mla    = (1/6)*sum(laXphi);     % Lame lambda(phi) integrated on reference element
    mmu    = (1/6)*sum(muXphi);     % Lame mu(phi) integrated on reference element
    
    aa(k,1,1,:) = mmu*(2*slocxx(:)*2+  slocyy(:)*2) + mla*(slocxx(:)*2);
    aa(k,2,2,:) = mmu*(  slocxx(:)*2+2*slocyy(:)*2) + mla*(slocyy(:)*2);
    aa(k,1,2,:) = mmu*(  slocyx(:)*2              ) + mla*(slocxy(:)*2);
    aa(k,2,1,:) = mmu*(  slocxy(:)*2              ) + mla*(slocyx(:)*2);
    
    % rhs
    muXphiXepsXphi = (muXphi')*(eps0*phiMatrix0ele(1,:));
    laXphiXepsXphi = (laXphi')*((eps0+eps0)*phiMatrix0ele(1,:));
    
    RHScoef = (1/24)*(2*(sum(diag(muXphiXepsXphi))+sum(muXphiXepsXphi(:))) ...
        +sum(diag(laXphiXepsXphi))+sum(laXphiXepsXphi(:)));
    
    sigmaXphi = sigma0*phiMatrix0ele(4,:);
    
    bElav1( k,1,:) = RHScoef* 2*slocx-sigmaXphi*clocx ;
    bElav1( k,2,:) = RHScoef* 2*slocy-sigmaXphi*clocy ;
    
    bElavF(k,1,:) = 1/6*mean(fMod(e2p(k,:)))*edet;
    bElavF(k,2,:) = 1/6*mean(fMod(npoint+e2p(k,:)))*edet;
end
%% Neumann boundary
neumann = [];
for k=1:nelement
    nodesDir = idp(e2p(k,:)) == 1;
    nodesNeu = idp(e2p(k,:)) == 2;
    if any(nodesNeu) && (sum(nodesDir + nodesNeu) >= 2)
        neumannNew = [e2p(k, nodesDir) e2p(k, nodesNeu)];
        neumann = [neumann; neumannNew];
    end
end





% % TODO only works for this special square example with x==-2 and x==2
bEla2 = zeros(2*npoint, 1);
bEla3 = zeros(2*npoint, 1);
for k=1:size(neumann,1)
    valueLength = abs(y(neumann(k,1)) - y(neumann(k,2))) / 2;
    value2 = mean(phi(neumann(k,:),4));
    value3_1 = mean((lambda+mu)*phi((neumann(k,1)),:)');
    value3_2 = mean(phi(neumann(k,:),1));
    if x(neumann(k,1)) == 2
        sign = 1;
    elseif x(neumann(k,1)) == -2
        sign = -1;
    else
        error('Not implemented :(');
    end
    bEla2(neumann(k,:)) = bEla2(neumann(k,:)) + sign*sigma0*value2*valueLength;
    bEla3(neumann(k,:)) = bEla3(neumann(k,:)) - sign*2*eps0*value3_1*value3_2*valueLength;
end







% iin     = zeros(size(neumann,1),2,2);
% jjn     = zeros(size(neumann,1),2,2);
% bElav2  = zeros(size(neumann,1),2,2);
% bElav3  = zeros(size(neumann,1),2,2);
% for k=1:size(neumann,1)
%     for l=1:2
%         iin( k,l,:) = (l-1)*npoint + neumann(k,:);
%         jjn( k,l,:) = ones(1,2);
%     end
%
%     val1_phi1 = phi(neumann(k,:),1);
%     val1_phi4 = phi(neumann(k,:),4);
%     val2 = mean((lambda+mu)*phi((neumann(k,:)),:)');
%     val3 = abs(y(neumann(k,1)) - y(neumann(k,2)))*[2 1; 1 2]/6;
%     if x(neumann(k,1)) == 2
%         sign = 1;
%     elseif x(neumann(k,1)) == -2
%         sign = -1;
%     else
%         error('Not implemented :(');
%     end
%
%     bElav2(k,1,:) = sign*sigma0*val1_phi4'*val3;
%     bElav3(k,1,:) = -sign*2*eps0*val2*val1_phi1'*val3;
% end
% bEla2 = sparse(iin(:),jjn(:),bElav2(:));
% bEla2 = [bEla2; zeros(2*npoint-length(bEla2),1)];
% bEla3 = sparse(iin(:),jjn(:),bElav3(:));
% bEla3 = [bEla3; zeros(2*npoint-length(bEla3),1)];






% iin     = zeros(size(neumann,1),2,2);
% jjn     = zeros(size(neumann,1),2,2);
% bElav2  = zeros(size(neumann,1),2,2);
% bElav3  = zeros(size(neumann,1),2,2);
% for k=1:size(neumann,1)
%     for l=1:2
%         iin( k,l,:) = (l-1)*npoint + neumann(k,:);
%         jjn( k,l,:) = ones(1,2);
%     end
%
%     val1_phi1 = phi(neumann(k,:),1);
%     val1_phi4 = phi(neumann(k,:),4);
%     val2 = mean((lambda+mu)*phi((neumann(k,:)),:)');
%     val3 = abs(y(neumann(k,1)) - y(neumann(k,2)))*[2 1; 1 2]/6;
%     if x(neumann(k,1)) == 2
%         sign = 1;
%     elseif x(neumann(k,1)) == -2
%         sign = -1;
%     else
%         error('Not implemented :(');
%     end
%
%     bElav2(k,1,:) = sign*sigma0*val1_phi4'*val3;
%
%     mat3_1 = abs(y(neumann(k,1)) - y(neumann(k,2)))*[3 1; 1 1]/12;
%     mat3_2 = abs(y(neumann(k,1)) - y(neumann(k,2)))*[1 1; 1 3]/12;
%     for m=1:6
%         alpha1 = phi(neumann(k,:),1);
%         alpham = phi(neumann(k,:),m);
%         bElav3(k,1,1) = bElav3(k,1,1) - sign*2*eps0*(lambda(m)+mu(m)) * alpha1'*mat3_1*alpham;
%         bElav3(k,1,2) = bElav3(k,1,2) - sign*2*eps0*(lambda(m)+mu(m)) * alpha1'*mat3_2*alpham;
%     end
% end
% bEla2 = sparse(iin(:),jjn(:),bElav2(:));
% bEla2 = [bEla2; zeros(2*npoint-length(bEla2),1)];
% bEla3 = sparse(iin(:),jjn(:),bElav3(:));
% bEla3 = [bEla3; zeros(2*npoint-length(bEla3),1)];











A    = sparse(ii(:),jj(:),aa(:));      % system matrix
bEla1 = sparse(iiv(:), jjv(:), bElav1(:));
bElaF = sparse(iiv(:), jjv(:), bElavF(:));
bEla = bEla1 + bElaF + bEla2 + bEla3;

% compute solution on Omega\Dirichlet
id    = ~[idp(1:npoint)==1;idp(1:npoint)==1]; % find indices of Dirichlet boundary
u     = zeros(2*npoint,1);
u(id) = A(id,id)\bEla(id);

%% Solve elasticity 2

addpath('./PDE');

coordinates = [x y];
elements3 = e2p;
DirichletNodes = find(idp == 1);
% Assembly
ANew = sparse(2*size(coordinates,1),2*size(coordinates,1));
for j = 1:size(elements3,1)
    I = 2*elements3(j,[1,1,2,2,3,3]) -[1,0,1,0,1,0];
    ANew(I,I) = ANew(I,I) +stima3(coordinates(elements3(j,:),:),lambda(1),mu(1));
end
pert = [linspace(2,2*npoint,npoint) - 1, linspace(2,2*npoint,npoint)];
ANew = ANew(pert, pert);
% Volume forces
bNew = zeros(2*size(coordinates,1),1);
for j = 1:size(elements3,1)
    I = 2*elements3(j,[1,1,2,2,3,3]) -[1,0,1,0,1,0];
    fs = [mean(fMod(elements3(j,:),:)); mean(fMod(npoint+elements3(j,:),:))];
    %     fs = fEval(sum(coordinates(elements3(j,:),:))/3)';
    %     fs = [f1(sum(coordinates(elements3(j,:),:))/3); f2(sum(coordinates(elements3(j,:),:))/3)];
    %     fs = [mean([f1(coordinates(elements3(j,1),:)), f1(coordinates(elements3(j,2),:)), f1(coordinates(elements3(j,3),:))]); mean([f2(coordinates(elements3(j,1),:)), f2(coordinates(elements3(j,2),:)), f2(coordinates(elements3(j,3),:))])];
    bNew(I) = bNew(I) +det([1,1,1;coordinates(elements3(j,:),:)'])*[fs;fs;fs]/6;
end
bNew = bNew(pert);

% % Dirichlet conditions
% M = repmat(eye(2), length(DirichletNodes), 1);
% W = zeros(2*length(DirichletNodes), 1);
% B = sparse(size(W,1),2*size(coordinates,1));
% for k = 0:1
%     for l = 0:1
%         B(1+l:2:size(M,1),2*DirichletNodes-1+k) = diag(M(1+l:2:size(M,1),1+k));
%     end
% end
% mask = find(sum(abs(B)'));
% ANew = [ANew, B(mask,:)'; B(mask,:), sparse(length(mask),length(mask))];
% bNew = [bNew;W(mask,:)];
%
% % Calculating the solution
% uNew = ANew \ bNew;
% uNew = uNew(1:2*size(coordinates,1));

uNew = zeros(2*npoint, 1);
DirichletNodes2 = 2*repmat(DirichletNodes, 2, 1) - [zeros(length(DirichletNodes), 1); ones(length(DirichletNodes), 1)];
FreeNodes2 = setdiff(1:(2*npoint), DirichletNodes2);
uNew(FreeNodes2) = ANew(FreeNodes2,FreeNodes2) \ bNew(FreeNodes2);

uNew = uNew(pert);

%% Check with the solution which should be optimal


gridDelaunay = delaunay(x,y);
fig1 = figure();
trimesh(gridDelaunay,x,y,uExpected(1:npoint));
fig2 = figure();
trimesh(gridDelaunay,x,y,u(1:npoint));
zlim([-0.2,1])

uDiff = uExpected - u;
sqrt(uDiff'*matrices.H1scal2D*uDiff)

fig3 = figure();
trimesh(gridDelaunay,x,y,uDiff(1:npoint));


% SaveTightFigure(fig1, 'Dir1');
% SaveTightFigure(fig2, 'Dir2');
% SaveTightFigure(fig3, 'Dir3');

SaveTightFigure(fig1, 'Neu1');
SaveTightFigure(fig2, 'Neu2');
SaveTightFigure(fig3, 'Neu3');

