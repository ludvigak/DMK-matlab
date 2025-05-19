clear

set(groot, 'defaultTextInterpreter','latex');

N = 5000;
[filepath,~,~] = fileparts(mfilename('fullpath'));
datafile = sprintf("%s/data/parameter_data_N%d.mat", filepath, N);

ds = load(datafile, 'data_dict');
data_dict = ds.data_dict;

% Make tolerance table
make_tol_table([kernels.rotlet_pswf.name, kernels.rotlet_ewald.name], data_dict)
make_tol_table([kernels.stokeslet_pswf.name, kernels.stokeslet_hasimoto.name], data_dict)
make_tol_table([kernels.stresslet_pswf.name, kernels.stresslet_hasimoto.name], data_dict)

% Make c-p plots
saveplots = true;
periodic = false;

make_cp_plots(kernels.rotlet_pswf.name, periodic, data_dict, saveplots);
make_cp_plots(kernels.stokeslet_pswf.name, periodic, data_dict, saveplots);
make_cp_plots(kernels.stresslet_pswf.name, periodic, data_dict, saveplots);

function make_tol_table(name_doublet, data_dict)
    tol_data = [];
    tol_list = 10.^(-3:-3:-12);
    for tol=tol_list
        tol_row = struct();
        tol_row.tol = tol;
        for name_idx=[1, 2]
            name = name_doublet(name_idx);
            input_data = data_dict{name};
            input_data.relerr_l2 = input_data.err_l2 ./ input_data.u_l2;
            for periodic = [false true]
                data =  input_data(input_data.periodic==periodic, :);
                if contains(name, 'PSWF')
                    shape = 'c';
                    suffix = '_pswf';
                else
                    shape = 'sigma';
                    suffix = '_gauss';
                end
                if periodic
                    suffix = [suffix '_per'];
                else
                    suffix = [suffix '_free'];
                end
                % Find best params to get below a given tolerance
                mask = data.relerr_l2 < tol;
                if ~any(mask)
                    subdata = data(1,:);
                    subdata.p = NaN;
                    subdata.(shape) = NaN;
                else
                    subdata = data(mask, :);
                end
                best = sortrows(subdata, {'p', shape, 'relerr_l2'});
                best = best(1,:);
                tol_row.([shape suffix]) = best.(shape);
                tol_row.(['p' suffix]) = best.p;
                %tol_row.(['err' suffix]) = best.relerr_l2;
            end
        end
        if isempty(tol_data)
            tol_data = tol_row;
        else
            tol_data(end+1) = tol_row;
        end
    end
    tol_data = struct2table(tol_data);

    % Collect max values over periodic and free/space
    max_data = table();
    max_data.tol = tol_data.tol;
    max_data.c_pswf = max(tol_data.c_pswf_free, tol_data.c_pswf_per);
    max_data.p_pswf = max(tol_data.p_pswf_free, tol_data.p_pswf_per);
    max_data.N1_pswf = 2*(max_data.c_pswf * 3/pi - 1) + 1;
    max_data.Nper_pswf = 2*(ceil(max_data.c_pswf / (2*pi)) - 1) + 1;

    max_data.sigma_gauss = max(tol_data.sigma_gauss_free, tol_data.sigma_gauss_per);
    max_data.p_gauss = max(tol_data.p_gauss_free, tol_data.p_gauss_per);
    Kmax = 2./max_data.sigma_gauss.^2;
    nf = 2*Kmax / (2*pi/3);
    max_data.N1_gauss = 2*(nf - 1) + 1;
    max_data.Nper_gauss = 2*(ceil(Kmax / (2*pi)) - 1) + 1;

    disp('----------')
    disp(name_doublet)
    disp(max_data)

    
    % Print Latex table
    header = '$\epsilon$ & $\frac{3}{\pi}c$ & $p$ & $N_1$ & $N_{\rm{per}}$ & $\frac{6}{\pi}\sigma^{-2}$ & $p^{\rm G}$ &  $N_1^{\rm G}$ & $N_{\rm{per}}^{\rm G}$ \\';
    body = sprintf('$10^{% 3d}$ & % 3d & % 3d & % 3d & % 3d & % 3.0f & % 3d & % 4.0f & % 3d \\\\\n', ...
                   [log10(max_data.tol) max_data.c_pswf*3/pi max_data.p_pswf max_data.N1_pswf max_data.Nper_pswf 6./(pi*max_data.sigma_gauss.^2) max_data.p_gauss max_data.N1_gauss max_data.Nper_gauss]');
    tbegin = '\begin{tabularx}{\textwidth}{Y|YYYY|YYYY}';
    tend = '\end{tabularx}';
    tabular = sprintf('%s\n\\hline\n%s\n\\hline\n%s\n\\hline\n%s\n', tbegin, header, body, tend)
    kname = split(name, ' ');
    kname = kname(1);
    filename = sprintf('tab/%s_table.tex', kname);
    fh = fopen(filename, 'w');
    fprintf(fh, '%s', tabular);
    fclose(fh);
    disp(['Wrote ' filename])
end

function make_cp_plots(name, periodic, data_dict, saveplots)
    data = data_dict{name};
    data =  data(data.periodic==periodic, :);
    data.relerr_l2 = data.err_l2 ./ data.u_l2;
    if contains(name, 'PSWF')
        shape = 'c';
    else
        shape = 'sigma';
    end

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

    % ==== PLOTS
    error_label = 'Relative $l_2$ error $E$';
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
    Ppc = polyfit(x(fit_mask), y(fit_mask), 1)
    %axis([0 inf 0 inf])
    grid on
    hold on
    xx = linspace(0, max(best_p_data.(shape))+5);
    plot(xx, polyval(Ppc, xx))
    xlabel(xlab)
    ylabel('p')
    axis equal
    axis([0 60 0 60])
    % P(2) always negative
    legend('Data', ...
           sprintf('$p = %.2f c %.2f$', Ppc(1), Ppc(2)), location='SouthEast',...
          interpreter='latex')
    setup_fig();

    sfigure(4);
    clf
    z = best_p_data.relerr_l2;
    semilogy(x, z, '.k')
    Pc = polyfit(x(fit_mask), log10(z(fit_mask)), 1)
    hold on
    semilogy(xx, 10.^( polyval(Pc, xx) ), '-')
    grid on
    xlim([0 45])
    xlabel(xlab)
    ylabel(error_label)
    % P(2) always negative
    if Pc(2) > 0
        plusstr = '+';
    else
        plusstr = '';
    end
    legend('Data', ...
           sprintf('$\\log_{10} E = %.2f c %s %.2f $', Pc(1), plusstr, Pc(2)),...
          interpreter='latex')

    sfigure(5);
    clf
    semilogy(best_p_data.p, best_p_data.relerr_l2, '.k')
    p_mask = (best_p_data.relerr_l2 < 1e-3) & (best_p_data.relerr_l2 > 1e-13);
    Pp = polyfit(best_p_data.p(p_mask), log10(best_p_data.relerr_l2(p_mask)), 1);
    hold on
    pplot = linspace(0, 60);
    semilogy(pplot, 10.^polyval(Pp, pplot))
    xlabel('$p$')
    ylabel(error_label)
    xlim([0 inf])
    grid on
    if Pp(2) > 0
        plusstr = '+';
    else
        plusstr = '';
    end
    legend('Data', ...
           sprintf('$\\log_{10} E = %.2f p %s %.2f $', Pp(1), plusstr, Pp(2)),...
          interpreter='latex')

    %

    % Equal y-scale for all error figs
    for f=[1 2 4 5]
        sfigure(f);
        ylim([1e-16 1])
        a = gca();
        a.YAxis.TickValues = 10.^(-15:5:0);
        setup_fig();
    end

    filename = strrep(lower(name), ' ', '_');
    if periodic
        filename = [filename '_periodic'];
    end
    if saveplots
        write_fig(1, sprintf('fig/%s_c_errors' , filename))
        write_fig(2, sprintf('fig/%s_p_errors' , filename))
        write_fig(3, sprintf('fig/%s_c_vs_p'   , filename))
        write_fig(4, sprintf('fig/%s_c_err_fit', filename))
        write_fig(5, sprintf('fig/%s_p_err_fit', filename))
    end
end

