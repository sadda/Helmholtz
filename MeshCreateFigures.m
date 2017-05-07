function figures = MeshCreateFigures(TriInfo)
    % Plots figures for mesh creation
    
    x        = TriInfo.x;
    y        = TriInfo.y;
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    npoint   = TriInfo.npoint;
    sizePhi  = TriInfo.sizePhi;
    
    % Plot random phi (with prescribed values)
    phi = rand(npoint,sizePhi);
    phi = phi ./ repmat(sum(phi,2), 1, sizePhi);
    phi = phi(TriInfo.phiRowsFree,:);
    phi = ProlongPhi(phi, TriInfo);
    gridDelaunay = delaunay(x,y);
    fig1 = figure();
    trimesh(gridDelaunay,x,y,phi(:,end));
    
    % Plot prescribed (not free) nodes
    fig2 = figure();
    hold off;
    for k=1:nelement
        patch(x(e2p(k,1:3))',y(e2p(k,1:3))','w');
        hold on;
    end
    for k=1:sizePhi
        scatter(x(phi(:,k)==1), y(phi(:,k)==1), 100, 'filled');
    end
    
    % Plot free nodes
    fig3 = figure();
    hold off;
    for k=1:nelement
        patch(x(e2p(k,1:3))',y(e2p(k,1:3))','w');
        hold on;
    end
    scatter(x(TriInfo.phiRowsFree), y(TriInfo.phiRowsFree), 100, 'filled');
    
    % Collect the figures into a struct
    figures = {fig1, fig2, fig3};
end


