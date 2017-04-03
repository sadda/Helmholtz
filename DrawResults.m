function DrawResults(phi, Theta, TriInfo, picName, data)
    
    phiProlonged = ProlongPhi(phi, TriInfo);
    
    minX      = min(TriInfo.x);
    maxX      = max(TriInfo.x);
    minY      = min(TriInfo.y);
    maxY      = max(TriInfo.y);
    
    %% Draw phi
    
    fig = figure('visible', 'off');
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiProlonged*(1:TriInfo.sizePhi)');
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    % title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', data.totalArea, data.J1, data.J1 / data.totalArea));
    view(2);
    shading interp;
    colormap(colormap(gray));
    % colormap jet;
    caxis([1 4]);
    c = colorbar;
    c.Ticks = [1 2 3 4];
    c.TickLabels = {'Ge', 'SiN', 'SiO_2', 'air'};
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    saveas(fig, fullfile(picName, sprintf('Phi_Ref%d_%d.jpg', data.meshIndex-1, data.iterIn)));
    
    %% Draw Theta
    
    fig = figure('visible', 'off');
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, Theta);
    % title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', data.totalArea, data.J1, data.J1 / data.totalArea));
    view(2);
    shading interp;
    colormap jet;
    colorbar;
    hold on;
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    if data.drawFixGe
        % Draw boundary edges for fixGe
        xEle     = TriInfo.x(TriInfo.e2p);
        yEle     = TriInfo.y(TriInfo.e2p);
        fixGeEle = sum(data.fixGe(TriInfo.e2p), 2) == 3;
        xDraw = xEle(fixGeEle,:);
        xDraw = xDraw(:,[1 1 2 2 3 3]);
        xDraw = reshape(xDraw, [], 2);
        yDraw = yEle(fixGeEle,:);
        yDraw = yDraw(:,[1 1 2 2 3 3]);
        yDraw = reshape(yDraw, [], 2);
        sortIndex          = (xDraw(:,1) > xDraw(:,2)) | (xDraw(:,1) == xDraw(:,2) & yDraw(:,1) > yDraw(:,2));
        xDraw(sortIndex,:) = xDraw(sortIndex,[2 1]);
        yDraw(sortIndex,:) = yDraw(sortIndex,[2 1]);
        xyDraw   = sortrows([xDraw yDraw]);
        xyDiff1  = diff(xyDraw);
        xyDiff2  = [ones(1, 4); xyDiff1];
        xyDiff3  = [xyDiff1; ones(1, 4)];
        xyUnique = sum(abs(xyDiff2),2)~=0 & sum(abs(xyDiff3),2)~=0;
        if ~isempty(xyDraw)
            xyDraw   = xyDraw(xyUnique,:);
            plot3(xyDraw(:,1:2)', xyDraw(:,3:4)', 2*max(Theta)*ones(size(xyDraw,1),2)', 'k--');
        else
            scatter3(TriInfo.x(data.fixGe), TriInfo.y(data.fixGe), 2*max(Theta)*ones(sum(data.fixGe),1), 'k', 'filled');
        end
        % Draw boundary edges for fixGePrev
        xEle     = TriInfo.x(TriInfo.e2p);
        yEle     = TriInfo.y(TriInfo.e2p);
        fixGeEle = sum(data.fixGePrev(TriInfo.e2p), 2) == 3;
        xDraw = xEle(fixGeEle,:);
        xDraw = xDraw(:,[1 1 2 2 3 3]);
        xDraw = reshape(xDraw, [], 2);
        yDraw = yEle(fixGeEle,:);
        yDraw = yDraw(:,[1 1 2 2 3 3]);
        yDraw = reshape(yDraw, [], 2);
        sortIndex          = (xDraw(:,1) > xDraw(:,2)) | (xDraw(:,1) == xDraw(:,2) & yDraw(:,1) > yDraw(:,2));
        xDraw(sortIndex,:) = xDraw(sortIndex,[2 1]);
        yDraw(sortIndex,:) = yDraw(sortIndex,[2 1]);
        xyDraw   = sortrows([xDraw yDraw]);
        xyDiff1  = diff(xyDraw);
        xyDiff2  = [ones(1, 4); xyDiff1];
        xyDiff3  = [xyDiff1; ones(1, 4)];
        xyUnique = sum(abs(xyDiff2),2)~=0 & sum(abs(xyDiff3),2)~=0;
        if ~isempty(xyDraw)
            xyDraw   = xyDraw(xyUnique,:);
            plot3(xyDraw(:,1:2)', xyDraw(:,3:4)', 2*max(Theta)*ones(size(xyDraw,1),2)', 'k');
        else
            scatter3(TriInfo.x(data.fixGePrev), TriInfo.y(data.fixGePrev), 2*max(Theta)*ones(sum(data.fixGePrev),1), 'k');
        end
    end
    
    saveas(fig, fullfile(picName, sprintf('Theta_Ref%d_%d.jpg', data.meshIndex-1, data.iterIn)));
end