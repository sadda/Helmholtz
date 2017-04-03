function fig = PlotFunction(f, TriInfo, data)
    
    isPhi = 0;
    if nargin < 3
        data = [];
    end
    if size(f,1) > TriInfo.npoint
        f     = reshape(f, [], TriInfo.sizePhi);
        isPhi = 1;
    end
    % Reduced phi
    if size(f,1) == TriInfo.npointRed && size(f,2) == TriInfo.sizePhi
        f     = ProlongPhi(f, TriInfo);
        isPhi = 1;
    end
    % Theta
    if size(f,1) == sum(~TriInfo.idp)
        fAux = zeros(TriInfo.npoint, size(f,2));
        fAux(~TriInfo.idp,:) = f;
        f = fAux;
    end
    if size(f,2) == TriInfo.sizePhi
        f     = f*(1:TriInfo.sizePhi)';
        isPhi = 1;
    end
    
    minX      = min(TriInfo.x);
    maxX      = max(TriInfo.x);
    minY      = min(TriInfo.y);
    maxY      = max(TriInfo.y);
    
    fig = figure();
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, f);
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    view(2);
    shading interp;
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    
    if isPhi
        colormap(colormap(gray));
        caxis([1 4]);
        c = colorbar;
        c.Ticks = [1 2 3 4];
        c.TickLabels = {'Ge', 'SiN', 'SiO_2', 'air'};
    else
        colormap(flipud(colormap(hot)));
        colorbar;
    end
        
    
    %     title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', totalArea, data.J1, data.J1 / totalArea));
end