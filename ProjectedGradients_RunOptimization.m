function [phi, t] = ProjectedGradients_RunOptimization(alpha, epsilon, TriInfo, Transformation, matrices, dirName, IterMax, drawResults, phi, tInitial)
    %RunOptimization
    %
    % In this main file the phase-field problem is solved using a
    % projected-gradient method, where the feasible set on which we project is
    % the Gibbs-simplex.
    
    if nargin < 10
        t = 0.5; % will be multiplied by 2
    else
        t = 0.5*tInitial;
    end
    
    if exist(dirName, 'dir')
        rmdir(dirName, 's');
    end
    mkdir(dirName);
    
    %% Set and save parameters
    [constants, material] = ObtainData(epsilon, alpha);
    sigma   = 1e-4;   % for Armijo line search
    tMin    = 1e-10;  % minimal step size
    tMax    = Inf;
    TOLabs  = 1e-5;
    TOLrel  = 1e-20;
    x       = TriInfo.x;
    y       = TriInfo.y;
    npoint  = TriInfo.npoint;
    e2p     = TriInfo.e2p;
    sizePhi = TriInfo.sizePhi;
    
    dataFile = 'Parameters.txt';
    fID_Data = fopen(dataFile,'w');
    fprintf(fID_Data,'epsilon %d\n alpha %d\n', constants.epsilon, constants.alpha);
    fprintf(fID_Data,'s %d\n l %d\n', constants.s, constants.l);
    fprintf(fID_Data,'sigma %d\n', sigma);
    fprintf(fID_Data,'tMin %d \n tMax %d\n', tMin, tMax);
    fprintf(fID_Data,'TOLrel  %d \n TOLabs = %d', TOLrel, TOLabs);
    fclose(fID_Data);
    
    %% Projected gradients
    iteration                      = 1;
    iterationData                  = 'Iteration_data.csv';
    res0                           = -Inf;
    res                            = Inf;
    resAll                         = nan(IterMax, 1);
    res2All                        = nan(IterMax, 1);
    lambda                         = zeros(size(phi));
    phiProj                        = phi;
    [JProj,~,JProj1,JProj2,JProj3] = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,0);
    
    while res > (1+res0)*TOLrel && res > TOLabs && iteration < IterMax && abs(t) >= tMin
        tic;
        phi           = phiProj;
        J             = JProj;
        J1            = JProj1;
        J2            = JProj2;
        J3            = JProj3;
        [~,gradient]  = ComputeData(phi,TriInfo,Transformation,matrices,constants,material,1);
        rieszGradient = ComputeRieszGradient(gradient, TriInfo, matrices);
        rieszGradient = reshape(rieszGradient,[],sizePhi);
        if iteration == 1
            rieszGradientNorm = sqrt(ComputePhiNormSquare(rieszGradient, TriInfo, matrices));
            cOptimality       = 1/rieszGradientNorm;
        end
        % Determine the next iterate
        t = min(3*t, tMax);
        [phiProj,t,lambda,JProj,JProj1,JProj2,JProj3,u,Theta] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin);
        % Compute the optimality (the same as in the loop with t=cOptimality)
        phiCheckNew                   = phi - cOptimality*rieszGradient;
        [phiCheck,~,~,iterationGibbs] = projection2Gibbs(phiCheckNew,phiProj,matrices,lambda,TriInfo);
        phiDiff                       = phi - phiCheck;
        res                           = sqrt(ComputePhiNormSquare(phiDiff, TriInfo, matrices));
        res2                          = sqrt(ComputePhiNormSquare(phi-phiProj, TriInfo, matrices));
        if iteration == 1
            res0 = res;
        end
        
        %% Print results and save iterations
        if iteration == 1 || mod(iteration, 20) == 0
            fprintf('%10s | %10s | %10s | %10s | %10s | %10s |\n', 'Iteration', 'StepSize', 'Objective', 'Residual', 'Residual2', 'Time');
        end
        elapsed = toc;
        fprintf('%10d | %3.4e | %3.4e | %3.4e | %3.4e | %3.4e |\n', iteration, t, J, res, res2, elapsed);
        
        if iteration == 1
            fID_Newton_p = fopen(iterationData,'w');
            fprintf(fID_Newton_p,'k t J res re2 iterationGibbs time J1 J2 J3');
            fprintf(fID_Newton_p,'\n %d %d %d %d %d %d %d %d %d %d %d %d %d',iteration-1,NaN,J,NaN,NaN,NaN,NaN,J1,J2,J3);
            fclose(fID_Newton_p);
        end
        
        fID_Newton_p = fopen(iterationData,'a+');
        fprintf(fID_Newton_p,'\n %d %d %d %d %d %d %d %d %d %d %d %d %d',iteration,t,J,res,res2,iterationGibbs,elapsed,J1,J2,J3);
        fclose(fID_Newton_p);
        
        filename = ['Iterate', num2str(iteration), '.mat'];
        save(filename, 'phi', 'u', 'Theta');
        movefile(filename,dirName);
        
        resAll(iteration)  = res;
        res2All(iteration) = res2;
        
        iteration = iteration+1;
    end
    
    fprintf('\n%13s %13s %13s %13s |\n', 'Objective', 'Objective1', 'Objective2', 'Objective3');
    fprintf('   %4.4e    %4.4e    %4.4e    %4.4e |\n', JProj, JProj1, JProj2, JProj3);
    fprintf('\n###########################################################################################################\n\n');
    
    save('phi','phi', 'u', 'Theta');
    save('DataAll');
    
    movefile('phi.mat',dirName);
    movefile('DataAll.mat',dirName);
    movefile(dataFile,dirName);
    movefile(iterationData,dirName);
    
    %% Draw all figures (after finishing optimization)
    if drawResults
        fileNumberSpace = 200;
    else
        fileNumberSpace = Inf;
    end
    % Load phi from the files
    if exist(strcat(dirName, '/IteratesAll.mat'), 'file') == 2
        load(strcat(dirName, '/IteratesAll.mat'));
    else
        fileNumber = 1;
        while exist(strcat(dirName, '/Iterate', int2str(fileNumber), '.mat'), 'file') == 2
            fileNumber = fileNumber + 1;
        end
        fileNumber = fileNumber - 1;
        load(strcat(dirName, '/Iterate', int2str(1), '.mat'));
        phiAll   = zeros(fileNumber, size(phi,1), size(phi,2));
        uAll     = zeros(fileNumber, size(u,1), size(u,2));
        ThetaAll = zeros(fileNumber, size(Theta,1), size(Theta,2));        
        for iteration=1:fileNumber
            load(strcat(dirName, '/Iterate', int2str(iteration), '.mat'));
            delete(strcat(dirName, '/Iterate', int2str(iteration), '.mat'));
            phiAll(iteration+1,:,:)   = phi;
            uAll(iteration+1,:,:)     = u;
            ThetaAll(iteration+1,:,:) = Theta;            
        end
        save(strcat(dirName, '/IteratesAll.mat'), 'phiAll', 'uAll', 'ThetaAll');
    end
    fileNumber = size(phiAll, 1);
    figureNumber = 1 + ceil(fileNumber / fileNumberSpace);
    if figureNumber <= 0 || figureNumber >= fileNumber
        figureNumber = fileNumber;
    end
    % Start plotting
    if figureNumber > 1 || ~exist('minIndex', 'var') % A bit hack here
        iterationAll = floor(linspace(1, fileNumber-1, figureNumber));
    else
        iterationAll = minIndex;
    end
    colormap jet;
    for iteration=iterationAll
        phi          = squeeze(phiAll(iteration+1,:,:));  
        u            = squeeze(uAll(iteration+1,:,:));
        Theta        = squeeze(ThetaAll(iteration+1,:,:));
        u            = u(:);
        Theta        = Theta(:);
        phiProlonged = ProlongPhi(phi, TriInfo);
        
        set(gcf,'Visible','off');
        filename = ['iterate', num2str(iteration), '.jpg'];
        clf;
        trisurf(e2p, x, y, phiProlonged*(1:sizePhi)');
        view(2);
        shading interp
        saveas(gcf, filename, 'jpg');
        movefile(filename,dirName);
        
        for i=1:sizePhi
            set(gcf,'Visible','off');
            filename = ['Phi', int2str(i), '_iterate', num2str(iteration),'.jpg'];
            clf;
            trisurf(e2p, x, y, phiProlonged(:,i));
            saveas(gcf, filename, 'jpg');
            movefile(filename, dirName);
        end
        
        set(gcf,'Visible','off');
        filename = ['Ux_iterate', num2str(iteration), '.jpg'];
        clf;
        trisurf(e2p, x, y, u(1:npoint));
        saveas(gcf, filename, 'jpg');
        movefile(filename,dirName);
        
        set(gcf,'Visible','off');
        filename = ['Uy_iterate', num2str(iteration),'.jpg'];
        clf;
        trisurf(e2p, x, y, u(npoint+1:end));
        saveas(gcf,filename,'jpg');
        movefile(filename,dirName);
        
        v      = matrices.Mloc2D\(matrices.Tr2D*u);
        vx     = v(1:npoint);
        vy     = v(npoint+1:2*npoint);
        tr_eps = (vx + vy)/2;
        set(gcf,'Visible','off');
        filename = ['U_BiaxialStrain', num2str(iteration), '.jpg'];
        clf;
        trisurf(e2p, x, y, tr_eps);
        saveas(gcf,filename,'jpg');
        movefile(filename,dirName);
        
        boundPhi = phiProlonged(:,1)>=0.3;
        boundX   = TriInfo.x(boundPhi);
        boundY   = TriInfo.y(boundPhi);
        bound    = boundary(boundX, boundY);
        set(gcf,'Visible','off');
        filename = ['Mode', num2str(iteration), '.jpg'];
        trisurf(e2p, x, y, Theta);
        view(2);
        shading interp;
        colormap jet;
        colorbar;
        hold on;
        plot3(boundX(bound), boundY(bound), 2*ones(size(bound)), 'k');
        saveas(gcf,filename,'jpg');
        movefile(filename,dirName);
    end
end

function [phiProj,t,lambda,JProj,JProj1,JProj2,JProj3,u,Theta] = PerformLineSearch(phi,J,rieszGradient,t,lambda,TriInfo,Transformation,matrices,constants,material,sigma,tMin)
    phiProj = phi;
    while true
        phiNew                                       = phi-t*rieszGradient;
        [phiProj,lambda]                             = projection2Gibbs(phiNew,phiProj,matrices,lambda,TriInfo);
        [JProj,~,JProj1,JProj2,JProj3,~,~,~,u,Theta] = ComputeData(phiProj,TriInfo,Transformation,matrices,constants,material,0);
        phiDiff                                      = phi-phiProj;
        normPhiDiffSquare                            = ComputePhiNormSquare(phiDiff, TriInfo, matrices);
        if JProj-J <= -(sigma/t)*normPhiDiffSquare || t < tMin
            break;
        else
            t = 0.5*t;
        end
    end
end


function phiNorm = ComputePhiNormSquare(phi, TriInfo, matrices)
    % It is prolonged by zero outside of the effective domain
    phiProlonged = ProlongPhi(phi(:), TriInfo) - TriInfo.phiProlongationVector(:);
    % It should really be divided by 2 and not by sqrt(2). Just in case that you want to spend another two hours by thinking about it.
    phiNorm = phiProlonged'*matrices.H1scal*phiProlonged / 2;
end
