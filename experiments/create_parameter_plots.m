clear

set(groot, 'defaultTextInterpreter','latex');

N = 5000;
[filepath,~,~] = fileparts(mfilename('fullpath'));
datafile = sprintf("%s/data/parameter_data_N%d.mat", filepath, N);

ds = load(datafile, 'data_dict');
data_dict = ds.data_dict;

% Make c-p plots
saveplots = true;
periodic = false;

make_cp_plots(kernels.rotlet_pswf.name, periodic, data_dict, saveplots);
make_cp_plots(kernels.stokeslet_pswf.name, periodic, data_dict, saveplots);
make_cp_plots(kernels.stresslet_pswf.name, periodic, data_dict, saveplots);

function make_cp_plots(name, periodic, data_dict, saveplots)
    data = data_dict{name};
    data =  data(data.periodic==periodic, :);
    data.relerr_l2 = data.err_l2 ./ data.u_l2;
    if contains(name, 'PSWF')
        shape = 'c';
    else
        shape = 'sigma';
    end
    % Find best params to get below a given tolerance
    tol_data = table();
    tol_list = 10.^(-1:-1:-15);
    for tol=tol_list
        mask = data.relerr_l2 < tol;
        if ~any(mask)
            continue
        end
        subdata = data(mask, :);
        best = sortrows(subdata, {'p', shape, 'relerr_l2'});
        row = best(1, :);
        row.tol = [tol];
        tol_data(end+1, :) = row;
    end
    tol_data

    % Find best c for each given p
    p_unique = unique(data.p)';
    p_data = table();
    for p_val=p_unique
        sub = data(data.p==p_val, :);
        [~, minidx] = min(sub.relerr_l2);
        p_data(end+1, :) = sub(minidx,:);
    end
    % Filter out best c values when non-unique
    best_p_data = table();
    for shape_val = unique(p_data.(shape))'
        mask = p_data.(shape)==shape_val;
        rows = p_data(mask, :);
        rows = sortrows(rows, 'relerr_l2');
        best_p_data(end+1, :) = rows(1, :);
    end



    % Print latex table
    %fprintf('%1.0e & % 3d & \n', [tol_data.tol round(tol_data.c*3/pi)]')


    % ==== PLOTS
    error_label = 'Relative $l_2$ error';
    if shape=='sigma'
        rescale = @(x) 1./x.^2;
        xlab = '$1/\sigma^2$';
    else
        rescale = @(x) x;
        xlab = '$c$';
    end

    sfigure(1);
    clf
    %p_plot = unique(data.p)';
    p_plot = 60:-5:5;
    for p = p_plot
        sub = data(data.p==p, :);
        sub = sortrows(sub, shape);
        semilogy(rescale(sub.(shape)), sub.relerr_l2, '-', LineWidth=1, displayname=sprintf('p=%d', p))
        hold on
    end
    grid on
    xlim([0 inf])
    xlabel(shape)
    ylabel(error_label)
    title('$p = 5:5:60$', interpreter='latex')

    sfigure(2);
    clf
    if shape=='sigma'
        shape_plot = unique(data.(shape))';
    else
        shape_plot = ((50:-5:5)*pi/3);
    end
    for s = shape_plot
        idx = find(data.(shape)==s);
        assert(isempty(idx)==false)
        sub = data(idx, :);
        sub = sortrows(sub, 'p');
        semilogy(sub.p, sub.relerr_l2, LineWidth=1)
        hold on
    end
    grid on
    xlabel('$p$')
    ylabel(error_label)
    title('$c = (5:5:50) \pi/3$', interpreter='latex')
    axis([0 60 1e-15 1])

    sfigure(3);
    clf
    x = rescale(best_p_data.(shape));
    y = best_p_data.p;
    plot(x, y, '.k')
    %axis equal
    fit_mask = (best_p_data.relerr_l2 > 1e-13) & (best_p_data.relerr_l2 < 1e-3);
    P = polyfit(x(fit_mask), y(fit_mask), 1)
    %axis([0 inf 0 inf])
    grid on
    hold on
    xx = linspace(0, max(x));
    plot(xx, polyval(P, xx))
    xlabel(xlab)
    ylabel('p')
    axis equal
    axis([0 60 0 60])
    % P(2) always negative
    legend('Data', ...
           sprintf('$p = %.2f c %.2f$', P(1), P(2)), location='SouthEast')
    setup_fig();

    sfigure(4);
    clf
    z = best_p_data.relerr_l2;
    semilogy(x, z, '.k')
    P2 = polyfit(x(fit_mask), log10(z(fit_mask)), 1)
    hold on
    semilogy(xx, 10.^( polyval(P2, xx) ), '-')
    grid on
    xlim([0 45])
    xlabel(xlab)
    ylabel(error_label)
    % P(2) always negative
    if P2(2) > 0
        plusstr = '+';
    else
        plusstr = '';
    end
    legend('Data', ...
           sprintf('$\\log_{10} p = %.2f c %s %.2f$', P2(1), plusstr, P2(2)))

    sfigure(5);
    clf
    semilogy(best_p_data.p, best_p_data.relerr_l2, '.k')
    xlabel('$p$')
    ylabel(error_label)
    xlim([0 inf])
    grid on

    % Equal y-scale for all error figs
    for f=[1 2 4 5]
        sfigure(f);
        ylim([1e-16 1])
        a = gca();
        a.YAxis.TickValues = 10.^(-15:5:0);
        setup_fig();
    end

    filename = strrep(lower(name), ' ', '_');
    if saveplots
        write_fig(1, sprintf('fig/%s_c_errors' , filename))
        write_fig(2, sprintf('fig/%s_p_errors' , filename))
        write_fig(3, sprintf('fig/%s_c_vs_p'   , filename))
        write_fig(4, sprintf('fig/%s_c_err_fit', filename))
        write_fig(5, sprintf('fig/%s_p_err_fit', filename))
    end
end

