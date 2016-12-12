function [phi, t, Theta, dataEigen] = ProjectedGradients2(TriInfo, Transformation, matrices, material, constants, dirName, IterMax, drawResults, phi, tInitial, dataEigen)
    
    if nargin < 11 || isempty(dataEigen)
        dataEigen = containers.Map('KeyType','double','ValueType','any');
    end
    dataEigen(-1) = 0;
    
    if nargin < 10
        t = 0.5;
    else
        t = 0.5*tInitial;
    end
    
    if exist(dirName, 'dir')
        rmdir(dirName, 's');
    end
    mkdir(dirName);
    
    %% Set parameters
    sigma   = 1e-4;   % for Armijo line search
    tMin    = 1e-10;  % minimal step size
    tMax    = 1e10;
    TOLabs  = 1e-5;
    x       = TriInfo.x;
    y       = TriInfo.y;
    npoint  = TriInfo.npoint;
    e2p     = TriInfo.e2p;
    sizePhi = TriInfo.sizePhi;
    
    %% Projected gradients
    iteration     = 1;
    iterationData = fullfile(dirName, 'IterationData.csv');
    res           = Inf;
    resAll        = nan(IterMax, 1);
    lambda        = zeros(size(phi));
    phiProj       = phi;
    [JProj,~,~,~,~,~,~,~,~,~,dataEigen] = ComputeData2(phi,TriInfo,Transformation,matrices,constants,material,0,dataEigen);
    
    while res > TOLabs && iteration < IterMax && abs(t) >= tMin
        tic;
        phi           = phiProj;
        J             = JProj;
        dataEigen(-1) = 0;
        [~,gradient,~,~,~,~,~,~,~,~,dataEigen]  = ComputeData2(phi,TriInfo,Transformation,matrices,constants,material,1,dataEigen);
        rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
        rieszGradient = reshape(rieszGradient,[],sizePhi);
        
        
        
        
        
        
        if SymmetryError(rieszGradient, TriInfo) >= 1e-8
            error('Symmetriation procedure failed');
        end
        rieszGradient = SymmetryCompute(rieszGradient, TriInfo);
        
        
        
        
        
        
        
        if iteration == 1
            rieszGradientNorm = sqrt(ComputePhiNormSquare(rieszGradient, TriInfo, matrices));
            cOptimality       = 1/rieszGradientNorm;
        end
        % Determine the next iterate
        t = min(2*t, tMax);
        [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen);
        % Compute the optimality (the same as in the loop with t=cOptimality)
        phiCheckNew                   = phi - cOptimality*rieszGradient;
        [phiCheck,~,~,iterationGibbs] = ProjectionGibbs(phiCheckNew,phiProj,matrices,lambda,TriInfo);
        phiDiff                       = phi - phiCheck;
        res                           = sqrt(ComputePhiNormSquare(phiDiff, TriInfo, matrices));
        
        %% Print results and save iterations
        if iteration == 1 || mod(iteration, 20) == 0
            fprintf('%10s | %10s | %10s | %10s | %10s |\n', 'Iteration', 'StepSize', 'Objective', 'Residual', 'Time');
        end
        elapsed = toc;
        fprintf('%10d | %3.4e | %3.4e | %3.4e | %3.4e |\n', iteration, t, J, res, elapsed);
        
        if iteration == 1
            fID_Newton_p = fopen(iterationData, 'w');
            fprintf(fID_Newton_p,'k t J res iterationGibbs time');
            fprintf(fID_Newton_p,'\n %d %d %d %d %d %d %d %d %d',iteration-1,NaN,J,NaN,NaN,NaN);
            fclose(fID_Newton_p);
        end
        
        fID_Newton_p = fopen(iterationData,'a+');
        fprintf(fID_Newton_p,'\n %d %d %d %d %d %d %d %d %d',iteration,t,J,res,iterationGibbs,elapsed);
        fclose(fID_Newton_p);
        
        filename = ['Phi_iterate', num2str(iteration), '.mat'];
        save(fullfile(dirName, filename), 'phi');
        
        resAll(iteration)  = res;
        iteration = iteration + 1;
    end
    
    [J, ~, J1, J2, J3, ~, ~, ~, ~, Theta] = ComputeData2(phi,TriInfo,Transformation,matrices,constants,material,0,dataEigen);
    
    fprintf('\n%13s %13s %13s %13s |\n', 'Objective', 'Objective1', 'Objective2', 'Objective3');
    fprintf('   %4.4e    %4.4e    %4.4e    %4.4e |\n', J, J1, J2, J3);
    fprintf('\n###########################################################################################################\n\n');
    
    save(fullfile(dirName, 'phi'), 'phi');
    save(fullfile(dirName, 'DataAll'));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Draw all figures (after finishing optimization)
    % Load phi from the files
    if exist(strcat(dirName, '/PhiAll.mat'), 'file') == 2
        load(strcat(dirName, '/PhiAll.mat'));
    else
        fileNumber = 1;
        while exist(strcat(dirName, '/Phi_iterate', int2str(fileNumber), '.mat'), 'file') == 2
            fileNumber = fileNumber + 1;
        end
        fileNumber = fileNumber - 1;
        phiNew = load(strcat(dirName, '/Phi_iterate', int2str(1), '.mat'));
        phiNew = phiNew.phi;
        phiAll = zeros(fileNumber, size(phiNew,1), size(phiNew,2));
        for iteration=1:fileNumber
            load(strcat(dirName, '/Phi_iterate', int2str(iteration), '.mat'));
            delete(strcat(dirName, '/Phi_iterate', int2str(iteration), '.mat'));
            phiAll(iteration+1,:,:) = phiNew;
        end
        save(strcat(dirName, '/PhiAll.mat'), 'phiAll');
    end
    if drawResults
        fileNumberSpace = 200;
        fileNumber = size(phiAll, 1);
        figureNumber = 1 + ceil(fileNumber / fileNumberSpace);
        if figureNumber <= 0 || figureNumber >= fileNumber
            figureNumber = fileNumber;
        end
        if figureNumber > 1 || ~exist('minIndex', 'var') % A bit hack here
            iterationAll = floor(linspace(1, fileNumber-1, figureNumber));
        else
            iterationAll = minIndex;
        end
        colormap jet;
        for iteration=iterationAll
            phiNew = squeeze(phiAll(iteration+1,:,:));
            phiProlonged = ProlongPhi(phiNew, TriInfo);
            
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['iterate', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, phiProlonged*(1:sizePhi)');
            view(2);
            shading interp
            saveas(gcf, filename, 'jpg');
            
            for i=1:sizePhi
                set(gcf,'Visible','off');
                filename = fullfile(dirName, ['Phi', int2str(i), '_iterate', num2str(iteration),'.jpg']);
                clf;
                trisurf(e2p, x, y, phiProlonged(:,i));
                saveas(gcf, filename, 'jpg');
            end
            
            [~, ~, ~, ~, ~, ~, ~, ~, u, Theta] = ComputeData2(phiNew,TriInfo,Transformation,matrices,constants,material,0,dataEigen);
            
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['Ux_iterate', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, u(1:npoint));
            saveas(gcf, filename, 'jpg');
            
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['Uy_iterate', num2str(iteration),'.jpg']);
            clf;
            trisurf(e2p, x, y, u(npoint+1:end));
            saveas(gcf,filename,'jpg');
            
            v      = matrices.Mloc2D\(matrices.Tr2D*u);
            vx     = v(1:npoint);
            vy     = v(npoint+1:2*npoint);
            tr_eps = (vx + vy)/2;
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['U_BiaxialStrain', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, tr_eps);
            saveas(gcf,filename,'jpg');
            
            boundPhi = phiProlonged(:,1)>=0.3;
            boundX   = TriInfo.x(boundPhi);
            boundY   = TriInfo.y(boundPhi);
            bound    = boundary(boundX, boundY);
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['Mode', num2str(iteration), '.jpg']);
            trisurf(e2p, x, y, Theta);
            view(2);
            shading interp;
            colormap jet;
            colorbar;
            hold on;
            plot3(boundX(bound), boundY(bound), 2*ones(size(bound)), 'k');
            saveas(gcf,filename,'jpg');
        end
    end
end

function [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen)
    phiProj = phi;
    while true
        phiNew                              = phi-t*rieszGradient;
        [phiProj,lambda]                    = ProjectionGibbs(phiNew,phiProj,matrices,lambda,TriInfo);
        dataEigen(-1)                       = t;
        [JProj,~,~,~,~,~,~,~,~,~,dataEigen] = ComputeData2(phiProj,TriInfo,Transformation,matrices,constants,material,0,dataEigen);
        phiDiff                             = phi-phiProj;
        normPhiDiffSquare                   = ComputePhiNormSquare(phiDiff, TriInfo, matrices);
        % You have to have isequal(phi, phiProj). The problem is the way into which you pass parameters into ComputeData2. If you compute Theta for phi,
        % the values in dataEigen will be updated. If you call the function again, it will make one more iteration for Theta computation, which will
        % result in a different Theta than before even though it was called with the same phi (but different dataEigen). Thus, the left-hand side
        % below may be extremely small but positive while the right-hand side may be zero.
        if JProj-J <= -(sigma/t)*normPhiDiffSquare || t < tMin || isequal(phi, phiProj)
            break;
        else
            t = 0.5*t;
        end
    end
    
    
    
    if SymmetryError(phiProj, TriInfo) >= 1e-8
        error('Symmetriation procedure failed');
    end
    phiProj = SymmetryCompute(phiProj, TriInfo);
    
    
end


function phiNorm = ComputePhiNormSquare(phi, TriInfo, matrices)
    % It is prolonged by zero outside of the effective domain
    phiProlonged = ProlongPhi(phi(:), TriInfo) - TriInfo.phiProlongationVector(:);
    % It should really be divided by 2 and not by sqrt(2). Just in case that you want to spend another two hours by thinking about it.
    phiNorm = phiProlonged'*matrices.H1scal*phiProlonged / 2;
end
