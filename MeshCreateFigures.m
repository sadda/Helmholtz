function figures = MeshCreateFigures(TriInfo)
    x        = TriInfo.x;
    y        = TriInfo.y;
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    npoint   = TriInfo.npoint;
    sizePhi  = TriInfo.sizePhi;
    % plot prescribed elements
    fig1 = figure();
    hold off;
    for k=1:nelement
        patch(x(e2p(k,1:3))',y(e2p(k,1:3))',TriInfo.ide(k));
        hold on;
    end
    hold on;
    scatter(x(TriInfo.idp==1), y(TriInfo.idp==1), 100, 'k', 'filled');
    hold off
    title(sprintf('npoint = %i, nelement = %i',npoint,nelement));
    fprintf('npoint = %i, nelement = %i \n',npoint,nelement);
    % plot random phi (with prescribed values)
    phi = rand(npoint,sizePhi);
    phi = phi ./ repmat(sum(phi,2), 1, sizePhi);
    phi = phi(TriInfo.phiRowsFree,:);
    phi = ProlongPhi(phi, TriInfo);
    gridDelaunay = delaunay(x,y);
    fig2 = figure();
    trimesh(gridDelaunay,x,y,phi(:,end));
    % plot prescribed nodes
    fig3 = figure();
    hold off;
    for k=1:nelement
        patch(x(e2p(k,1:3))',y(e2p(k,1:3))','w');
        hold on;
    end
    for k=1:sizePhi
        scatter(x(phi(:,k)==1), y(phi(:,k)==1), 100, 'filled');
    end
    % plot free nodes
    fig4 = figure();
    hold off;
    for k=1:nelement
        patch(x(e2p(k,1:3))',y(e2p(k,1:3))','w');
        hold on;
    end
    scatter(x(TriInfo.phiRowsFree), y(TriInfo.phiRowsFree), 100, 'filled');
    figures = {fig1, fig2, fig3, fig4};
end