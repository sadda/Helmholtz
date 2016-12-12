clear all;

plotResults = 0;
meshRefine  = 0;
load(['Results_Ref' int2str(meshRefine), '/DataAll.mat'], 'TriInfo', 'matrices');

npoint = TriInfo.npoint;
id     = ~TriInfo.idp(1:TriInfo.npoint)==1;

S     = matrices.GradSq(1:npoint,1:npoint);
M     = matrices.Mloc(1:npoint,1:npoint);
Sid   = S(id,id);
MSid1 = diag(1./sum(M(id,id),2))*S(id,id);
MSid2 = M(id,id) \ S(id,id);

MSid  = MSid2;

a  = max(TriInfo.x)-min(TriInfo.x);
b  = max(TriInfo.y)-min(TriInfo.y);
ab = a*b;

%% True values
maxNumber = 1;
maxIter   = 20;
sq1 = (1:maxNumber).^2;
sq2 = (1:maxNumber).^2;
sq1 = reshape(repmat(sq1, maxNumber, 1), [] , 1);
sq2 = reshape(repmat(sq2, 1, maxNumber), [] , 1);

lambdaTrue = sort((sq1+sq2)*pi^2/ab);
lambdaTrue = lambdaTrue(1:maxNumber);

%% Using eigs
[xEigs1, lambdaEigs1] = eigs(MSid1, maxNumber, 'SM');
lambdaEigs1 = sort(diag(lambdaEigs1));

tic;
[~, lambdaEigs2] = eigs(MSid2, maxNumber, 'SM');
timeEigs    = toc;

if plotResults
    [~, lambdaEigs0] = eigs(Sid, maxNumber, 'SM');
    
    lambdaEigs0 = sort(diag(lambdaEigs0));
    lambdaEigs2 = sort(diag(lambdaEigs2));

    fig = figure;
    set(groot,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|:|--|-.')
    plot(log10([lambdaTrue lambdaEigs0 lambdaEigs1 lambdaEigs2]));
    legend({'True', 'Eigs S', 'Eigs M^{-1}S, M lumped', 'Eigs M^{-1}S'});
    ylim([-4 4]);
    xlabel('n');
    ylabel('log10[lambda(n)]');
    saveas(fig, ['Results' int2str(meshRefine), '.jpg']);
end

%% Inverse iteration (more intelligent) with factorization (but without shift)

tic;
y      = rand(sum(id),1);
[L,U]  = lu(S(id,id));
for i=1:maxIter
    x      = y / norm(y);
    Mx     = M(id,id)*x;
    y      = U \ (L \ Mx);
    mu     = y'*x;
    lambda = 1/mu;
    res    = S(id,id)*x-lambda*Mx;
    if norm(res) <= 1e-8
        break
    end
end
iInverse1      = i;
lambdaInverse1 = lambda;
xInverse1      = x;
timeInverse1   = toc;

%% Inverse iteration with factorization and shift

tic;
shift = 0;
A      = MSid;
AOrig  = A;
A      = A - shift*eye(size(A));
y      = rand(size(A,1),1);
[L,U]  = lu(A);
for i=1:maxIter
    x      = y / norm(y);
    y      = U \ (L \ x);
    mu     = y'*x;
    lambda = shift + 1/mu;
    res    = AOrig*x-lambda*x;
    if norm(res) <= 1e-4
        break
    end
end
iInverse2      = i;
lambdaInverse2 = lambda;
xInverse2      = x;
timeInverse2   = toc;

%% Rayleigh iteration

tic;
A = MSid;
x = xInverse1;
for i=1:maxIter
    rho   = x'*A*x;
    z     = (A-rho*eye(size(A,1))) \ x;
    x     = z / norm(z);
    res   = A*x-rho*x;
    if norm(res) <= 1e-10
        break
    end
end
iRayleigh      = i;
lambdaRayleigh = rho;
xRayleigh      = x;
timeRayleigh   = toc;


%% Print data

fprintf('%10s | %10s | %10s | %10s\n', 'Method', 'Iterations', 'Residual', 'Time');
fprintf('%10s | %10s | %10s | %10.4f\n', 'Eigs', '-', '-', timeEigs);
fprintf('%10s | %10d | %10.3e | %10.4f\n', 'Inverse 1', iInverse1, abs(lambdaInverse1 - lambdaEigs2), timeInverse1);
fprintf('%10s | %10d | %10.3e | %10.4f\n', 'Inverse 2', iInverse2, abs(lambdaInverse2 - lambdaEigs2), timeInverse2);
fprintf('%10s | %10d | %10.3e | %10.4f\n', 'Rayleigh', iRayleigh, abs(lambdaRayleigh - lambdaEigs2), timeRayleigh);

