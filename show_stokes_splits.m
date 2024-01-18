% Show Fourier and real space decay of different Stokeslet splits
clear

clf
tol = 1e-9;
show_split(@kernels.stokeslet_hasimoto, tol);
show_split(@kernels.stokeslet_pswf, tol);
show_split(@kernels.stokeslet_pswf_sq, tol);
show_split(@kernels.stokeslet_exp4, tol);

for p=1:2
    subplot(1,2,p)
    xa = xlim();
    plot(xa, [tol tol], '--', displayname='tolerance')
    legend(location='southoutside')
end

function show_split(kernel_ref, tol)
    kernel = kernel_ref(tolerance=tol)
    level = 0;
    k = linspace(0, kernel.Kmax, 200);
    gamma_hat = kernel.fourier_scaling(k.^2, level);

    r = linspace(0.001, 1, 200);
    [Rdiag, Roffd] = kernel.real_decay(r, level);
    
    subplot(1,2,1)
    semilogy(k, abs(gamma_hat), '.-', displayname=kernel.name)
    hold on
    grid on
    ylim([-inf 1])
    xlabel('k')
    title('Fourier scaling')
            
    subplot(1,2,2)
    plot(r, abs(Roffd)+abs(Rdiag), displayname=kernel.name)
    hold on
    set(gca, yscale='log')
    grid on
    ylim([tol/1e2 1])
    xlabel('r')
    title('Real scaling')
end
