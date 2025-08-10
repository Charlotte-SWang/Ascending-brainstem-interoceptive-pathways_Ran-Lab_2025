function removeGrayMargin(f,xsize,ysize)
%     % remove gray margin of the figure;
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
    ax.Position = [0 0 1 1];
    f.Position = [0 0 xsize ysize];
    F = getframe(gcf);
    [X, Map] = frame2im(F);
end

