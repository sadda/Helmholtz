function [phi, t, Theta, dataEigen] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, IterMax, drawResults, phi, t, options, dataEigen)
    % Runs the projected gradient method as described in the enlosed text.
    
    if nargin < 12 || isempty(dataEigen) % Stores (selected) eigenvalues computed earlier. Used for shift for the eigenvalue computation.
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
    
    sigma   = 1e-4;              % For the Armijo line search
    tMin    = 1e-10;             % Minimal step size
    tMax    = 1e10;              % Maximal step size
    stopTol = 1e-7;              % Stopping tolerance
    
    %% Set the initial data
    
    x       = TriInfo.x;         % Coordinates x
    y       = TriInfo.y;         % Coordinates y
    npoint  = TriInfo.npoint;    % Number of nodes
    e2p     = TriInfo.e2p;       % Elements
    sizePhi = TriInfo.sizePhi;   % Number for phases
    
    iteration     = 1;
    iterationData = fullfile(dirName, 'IterationData.csv');
    res           = Inf;              % Current residual
    resAll        = nan(IterMax, 1);  % All residuals
    lambda        = zeros(size(phi)); % Multiplier for the projection onto the Gibbs simplex
    phiProj       = phi;
    
    options.computeG = 0;             % Determines whether the gradient will be computed (saves time if not). Is switched two times during each iteration
    [JProj,~,~,~,~,~,~,~,~,~,dataEigen] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
    
    while res > stopTol && iteration <= IterMax && abs(t) >= tMin
        %% Run the optimization
        
        tic;
        phi           = phiProj; % phi
        J             = JProj;   % Objective function
        dataEigen(-1) = 0;
        options.computeG = 1;
        
        % First compute gradient and then its Riesz representation
        [~,gradient,~,~,~,~,~,~,u,Theta,dataEigen]  = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
        rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
        rieszGradient = reshape(rieszGradient,[],sizePhi);
        if options.symmetrize
            rieszGradient = SymmetryCompute(rieszGradient, TriInfo, 1, 0, 1e-8);
        end
        
        % Determine the step size and iterate
        t = min(2*t, tMax);
        [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen,options);
        
        % Compute the optimality (with t=cOptimality)
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
    
    %% Recompute some data and save everything
    
    options.computeG   = 0;
    [J, ~, J1, J2, J3] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,options,dataEigen);
    
    fprintf('\n%13s %13s %13s %13s |\n', 'Objective', 'Objective1', 'Objective2', 'Objective3');
    fprintf('   %4.4e    %4.4e    %4.4e    %4.4e |\n', J, J1, J2, J3);
    fprintf('\n###########################################################################################################\n\n');
    
    save(fullfile(dirName, 'phi'), 'phi');
    save(fullfile(dirName, 'DataAll'));
    
    %% Load phi, u and Theta from individual files and merge them into one big file
    
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
    
    %% Draw all figures
    
    if drawResults
        % Select which iterations will be drawn
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
        
        for iteration=iterationAll
            phi          = squeeze(phiAll(iteration,:,:));
            u            = squeeze(uAll(iteration,:))';
            Theta        = squeeze(ThetaAll(iteration,:))';
            phiProlonged = ProlongPhi(phi, TriInfo);
            
            % Draw phi (all at once)
            fig = PlotFunction(phiProlonged, TriInfo, 0);
            filename = fullfile(dirName, ['PhiAll', num2str(iteration), '.jpg']);
            saveas(fig, filename, 'jpg');
            
            % Plot Theta
            fig = PlotFunction(Theta, TriInfo, 0);
            filename = fullfile(dirName, ['Theta', num2str(iteration), '.jpg']);
            saveas(fig, filename, 'jpg');
            
            colormap jet;
            % Draw phi (each phase separately)
            for i=1:sizePhi
                set(gcf,'Visible','off');
                filename = fullfile(dirName, ['Phi', int2str(i), '_iterate', num2str(iteration),'.jpg']);
                clf;
                trisurf(e2p, x, y, phiProlonged(:,i));
                saveas(gcf, filename, 'jpg');
            end
            
            % Draw ux
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['Ux_iterate', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, u(1:npoint));
            saveas(gcf, filename, 'jpg');
            
            % Draw uy
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['Uy_iterate', num2str(iteration),'.jpg']);
            clf;
            trisurf(e2p, x, y, u(npoint+1:end));
            saveas(gcf,filename,'jpg');
            
            % Draw biaxial strain
            v      = matrices.Mloc2D\(matrices.Tr2D*u);
            vx     = v(1:npoint);
            vy     = v(npoint+1:2*npoint);
            tr_eps = (vx + vy)/2;
            set(gcf,'Visible','off');
            filename = fullfile(dirName, ['U_BiaxialStrain', num2str(iteration), '.jpg']);
            clf;
            trisurf(e2p, x, y, tr_eps);
            saveas(gcf,filename,'jpg');
        end
    end
end

function [phiProj,t,lambda,JProj,dataEigen] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin,dataEigen,options)
    % Computes the step size based on the Armijo condition
    
    phiProj = phi;
    while true
        phiNew                              = phi-t*rieszGradient;
        [phiProj,lambda]                    = ProjectionGibbs(phiNew,phiProj,matrices,lambda,TriInfo);
        dataEigen(-1)                       = t;
        options.computeG                    = 0;
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
    % Computes the H^1 norm of phi.
    
    phiProlonged  = ProlongPhi(phi(:), TriInfo) - TriInfo.phiProlongationVector(:);
    phiNormSquare = phiProlonged'*matrices.H1scal*phiProlonged ;
end

function rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices)
    % Computes the Riesz gradient, thus pulls the gradient from (H^1)^* into H^1
    
    phiRowsFree6  = TriInfo.phiRowsFree6;
    rieszGradient = matrices.H1scal(phiRowsFree6,phiRowsFree6) \ gradient;
end



