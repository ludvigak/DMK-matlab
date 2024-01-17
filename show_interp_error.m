% Show kernel interpolation error for varying tolerances and p
clear
figure(1); clf

kernel_list = {@kernels.laplace_ewald, @kernels.stokeslet_hasimoto};
for kidx = 1:numel(kernel_list)
    kernel_ref = kernel_list{kidx};
    step = numel(kernel_list);
    figure(2 + step*(kidx-1))
    clf
    rng(1);
    plist = 2:65;
    tol_list = logspace(-1, -13, 100);
    tol_plot = 10.^[-3 -6 -9 -12];

    tol_list = fliplr(unique(sort([tol_list, tol_plot])));
    sources = rand(1, 100)-1/2;
    charges = rand(1, 100)-1/2;
    p_first = [];
    err_first = [];
    tol_achieved = [];
    for tol = tol_list
        kernel = kernel_ref(tolerance=tol);
        interp_err = estimate_interp_error(plist, kernel);
        first = find(interp_err < tol, 1);
        if ~isempty(first)
            p_first(end+1) = plist(first);
            tol_achieved(end+1) = tol;
            err_first(end+1) = interp_err(first);
        end
        if any(tol==tol_plot)
            semilogy(plist, interp_err, '.-', 'DisplayName', sprintf('tol=%.0e', tol))
            hold on
            grid on
            legend
            xlabel('p')
            ylabel('rel. error')
        end
    end

    for t=tol_plot
        plot(plist, t*plist.^0, '-k', displayname=sprintf('level=%.0e', t))
    end
    title([kernel.name " interpolation error"])

    figure(1)
    semilogx(tol_achieved, p_first, '.-', displayname=kernel.name)
    xlabel("Error")
    ylabel("Minimum p")
    hold on
    grid on
    legend
end
