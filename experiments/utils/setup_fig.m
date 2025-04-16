function setup_fig(args)
    arguments
        args.fig = 0
        args.width_cm = 7.3;
        args.height_cm = 5.8;
    end
    if args.fig
        f = sfigure(args.fig);
    else
        f = gcf();
    end
    % Setup figure for paper plot
    set(gca,'FontName','Times','FontSize',8)
    f.Units = 'centimeter';
    f.PaperUnits = 'centimeter';
    f.PaperPositionMode = 'auto';
    P = f.Position;
    f.Position = [P(1) P(2) args.width_cm args.height_cm];
end

