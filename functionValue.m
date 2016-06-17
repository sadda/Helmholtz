function [J, J1, J2, J3] = functionValue(phi,u,Theta,TriInfo,Transformation,matrices,constants)
    % FUNCTIONVALUE gives back the funcation value of the objective functional
    % J of the phase-field minimization problem
    %
    % Input: phi            - phase-field vector
    %        u              - state variable (dependent on phi)
    %        TriInfo        - Information of triangulation as nelement, npoint, etc
    %        matrices       - matrices of FEM
    %        constants      - values of constants in a structure (gamma, M, c)
    %
    % Output: J - function value at (phi,u)
    
    %% extract information of triangulation, material properties and constants
    Mloc    = matrices.Mloc;
    Id      = matrices.Id;
    
    GradSq = matrices.GradSq;
    
    alpha   = constants.alpha;
    epsilon = constants.epsilon;
    
    phi = ProlongPhi(phi, TriInfo);
    
    
    
    
    npoint                  = TriInfo.npoint;
    nelement                = TriInfo.nelement;
    e2p                     = TriInfo.e2p;
    id                      = ~TriInfo.idp(1:TriInfo.npoint)==1;
    sizePhi                 = TriInfo.sizePhi;
    constants.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
    constants.epsilonR(2:3) = [];
    epsilonR                = constants.epsilonR;
    
    %% Compute objective
    
    A = zeros(nelement,1,1,9);
    for k=1:nelement
        edet   = Transformation{k,1};
        slocx  = Transformation{k,9};
        slocy  = Transformation{k,10};
        mloc   = Transformation{k,11};
        
        ux   = u(e2p(k,:));
        uy   = u(e2p(k,:)+npoint);
        
        dxux = ux'*slocx';
        dyuy = uy'*slocy';
        
        A(k,1,1,:) = 2*(dxux+dyuy)/edet*mloc(:);
    end
    ii1D = TriInfo.indicesIPhi(:,1,1,:);
    jj2D = TriInfo.indicesJPhi(:,1,1,:);
    A    = sparse(ii1D(:), jj2D(:), A(:));
    
    J1  = -Theta'*A*Theta;
    
    
    
    
    
%     J1 = -(matrices.TrXD+matrices.TrYD)'*u;
    J2 = 0.5*(GradSq*phi(:))'*phi(:);
    J3 = 0.5*(Id'*phi(:) - (Mloc*phi(:))'*phi(:));
    
    J  = J1 + alpha*epsilon*J2 + alpha/epsilon*J3;
end