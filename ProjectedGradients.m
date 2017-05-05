function [phi, t, Theta, dataEigen, data] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, IterMax, drawResults, phi, t, options, dataEigen)
    
    if nargin < 12 || isempty(dataEigen)
        dataEigen = containers.Map('KeyType','double','ValueType','any');
    end
    dataEigen(-1) = 0;
    
    if exist(dirName, 'dir')
        rmdir(dirName, 's');
    end
    mkdir(dirName);
    u     = []; % Needed for parfor
    Theta = [];
    
    %% Set parameters
    sigma   = 1e-4;   % for Armijo line search
    tMin    = 1e-10;  % minimal step size
    tMax    = 1e10;
    TOLabs  = 1e-7;
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
    
    options.computeG = 0;
    [JProj,~,~,~,~,~,~,~,~,~,dataEigen] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
    
    while res > TOLabs && iteration <= IterMax && abs(t) >= tMin
        tic;
        phi           = phiProj;
        J             = JProj;
        dataEigen(-1) = 0;
        options.computeG = 1;
        [~,gradient,~,~,~,~,~,~,u,Theta,dataEigen]  = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
        rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
        rieszGradient = reshape(rieszGradient,[],sizePhi);
        if options.symmetrize
            rieszGradient = SymmetryCompute(rieszGradient, TriInfo, 1, 0, 1e-8);
        end
        % Determine the next iterate
        if options.method == 1
            t = min(3*t, tMax);
        else
            t = min(2*t, tMax);
        end
        [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen,options);
        % Compute the optimality (the same as in the loop with t=cOptimality)
        phiCheckNew                   = phi - constants.cOptimality*rieszGradient;
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
        filename = ['U_iterate', num2str(iteration), '.mat'];
        save(fullfile(dirName, filename), 'u');
        filename = ['Theta_iterate', num2str(iteration), '.mat'];
        save(fullfile(dirName, filename), 'Theta');
        
        resAll(iteration) = res;
        iteration         = iteration + 1;
    end
    
    options.computeG   = 0;
    [J, ~, J1, J2, J3] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
    
    data.J  = J;
    data.J1 = J1;
    
    fprintf('\n%13s %13s %13s %13s |\n', 'Objective', 'Objective1', 'Objective2', 'Objective3');
    fprintf('   %4.4e    %4.4e    %4.4e    %4.4e |\n', J, J1, J2, J3);
    fprintf('\n###########################################################################################################\n\n');
    
    save(fullfile(dirName, 'phi'), 'phi');
    save(fullfile(dirName, 'DataAll'));
    
    %% Draw all figures (after finishing optimization)
    % Load phi from the files
    if exist(strcat(dirName, '/PhiAll.mat'), 'file') == 2
        load(strcat(dirName, '/PhiAll.mat'));
        load(strcat(dirName, '/UAll.mat'));
        load(strcat(dirName, '/ThetaAll.mat'));
    else
        fileNumber = 1;
        while exist(strcat(dirName, '/Phi_iterate', int2str(fileNumber), '.mat'), 'file') == 2
            fileNumber = fileNumber + 1;
        end
        fileNumber = fileNumber - 1;
        
        phiAll   = zeros(fileNumber, size(phi,1), size(phi,2));
        uAll     = zeros(fileNumber, size(u,1), size(u,2));
        ThetaAll = zeros(fileNumber, size(Theta,1), size(Theta,2));
        for iteration=1:fileNumber
            load(strcat(dirName, '/Phi_iterate', int2str(iteration), '.mat'));
            delete(strcat(dirName, '/Phi_iterate', int2str(iteration), '.mat'));
            phiAll(iteration,:,:) = phi;
            load(strcat(dirName, '/U_iterate', int2str(iteration), '.mat'));
            delete(strcat(dirName, '/U_iterate', int2str(iteration), '.mat'));
            uAll(iteration,:,:) = u;
            load(strcat(dirName, '/Theta_iterate', int2str(iteration), '.mat'));
            delete(strcat(dirName, '/Theta_iterate', int2str(iteration), '.mat'));
            ThetaAll(iteration,:,:) = Theta;
        end
        save(strcat(dirName, '/PhiAll.mat'), 'phiAll');
        save(strcat(dirName, '/UAll.mat'), 'uAll');
        save(strcat(dirName, '/ThetaAll.mat'), 'ThetaAll');
    end
    if drawResults
        fileNumberSpace = 200;
        fileNumber      = size(phiAll, 1);
        figureNumber    = 1 + ceil(fileNumber / fileNumberSpace);
        if figureNumber <= 0 || figureNumber >= fileNumber
            figureNumber = fileNumber;
        end
        if figureNumber > 1 || ~exist('minIndex', 'var')
            iterationAll = floor(linspace(1, fileNumber, figureNumber));
        else
            iterationAll = minIndex;
        end
        colormap jet;
        for iteration=iterationAll
            phi          = squeeze(phiAll(iteration,:,:));
            u            = squeeze(uAll(iteration,:))';
            Theta        = squeeze(ThetaAll(iteration,:))';
            phiProlonged = ProlongPhi(phi, TriInfo);
            
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['iterate', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, phiProlonged*(1:sizePhi)');
            view(2);
            shading interp;
            saveas(gcf, filename, 'jpg');
            
            for i=1:sizePhi
                set(gcf,'Visible','off');
                filename = fullfile(dirName, ['Phi', int2str(i), '_iterate', num2str(iteration),'.jpg']);
                clf;
                trisurf(e2p, x, y, phiProlonged(:,i));
                saveas(gcf, filename, 'jpg');
            end
            
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

function [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen,options)
    
    phiProj = phi;
    while true
        phiNew                              = phi-t*rieszGradient;
        [phiProj,lambda]                    = ProjectionGibbs(phiNew,phiProj,matrices,lambda,TriInfo);
        dataEigen(-1)                       = t;
        options.computeG = 0;
        [JProj,~,~,~,~,~,~,~,~,~,dataEigen] = ComputeData(phiProj,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
        phiDiff                             = phi-phiProj;
        normPhiDiffSquare                   = ComputePhiNormSquare(phiDiff, TriInfo, matrices);
        if JProj-J <= -(sigma/t)*normPhiDiffSquare || t < tMin
            break;
        else
            t = 0.5*t;
        end
    end
    if options.symmetrize
        phiProj = SymmetryCompute(phiProj, TriInfo, 1, 0, 1e-8);
    end
end


function phiNormSquare = ComputePhiNormSquare(phi, TriInfo, matrices)
    
    % It is prolonged by zero outside of the effective domain
    phiProlonged  = ProlongPhi(phi(:), TriInfo) - TriInfo.phiProlongationVector(:);
    phiNormSquare = phiProlonged'*matrices.H1scal*phiProlonged ;
end
