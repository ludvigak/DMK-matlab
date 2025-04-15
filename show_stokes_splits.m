% Show Fourier and real space decay of different Stokeslet splits
clear

sfigure(1); clf
tol = 1e-6;
c_pswf=18;
show_split(kernels.stokeslet_hasimoto(tolerance=tol), tol, '.-');
show_split(kernels.stokeslet_pswf(c_pswf=c_pswf)    , tol, 'o-');
show_split(kernels.stokeslet_pswf_num(c_pswf=c_pswf)    , tol, '*-');
%show_split(kernels.stokeslet_pswf3(c_pswf=c_pswf)    , tol, 'p-');
% These are clearly inferior in real space:
%show_split(@kernels.stokeslet_pswf_sq, tol);
%show_split(@kernels.stokeslet_exp4, tol);

for p=1:2
    subplot(1,2,p)
    xa = xlim();
    plot(xa, [tol tol], '--', displayname='tolerance')
    legend(location='southoutside')
end

sfigure(2);clf
c_pswf_vec = 2:0.5:50;
plot_truncation(@kernels.stokeslet_pswf, c_pswf_vec);
plot_truncation(@kernels.stresslet_pswf, c_pswf_vec);

function plot_truncation(kernel_ref, c_pswf_vec)
    [Rtrunc, Ftrunc] = get_truncation(kernel_ref, c_pswf_vec);
    name = kernel_ref(tol=0.1).name;
    subplot(1, 2, 1)
    semilogy(c_pswf_vec, Rtrunc, displayname=name)
    ylim([1e-16 10])
    grid on
    hold on
    legend()
    subplot(1,2,2)
    semilogy(c_pswf_vec, Ftrunc, displayname=name)
    ylim([1e-16 10])
    grid on
    hold on
    legend()
end
function [Rtrunc, Ftrunc] = get_truncation(kernel_ref, c_pswf_vec)
    [Rtrunc, Ftrunc] = deal(zeros(size(c_pswf_vec)));
    for i=1:numel(c_pswf_vec)
        c = c_pswf_vec(i);
        kernel = kernel_ref(c_pswf=c);
        f = ones(1, kernel.dim_in);
        Rtrunc(i) = max(abs(kernel.reskernel([1 0 0], [0 0 0], f, 0)));
        [Rdiag, Roffd] = kernel.real_decay(1, 0);
        %Rtrunc(i) = max(abs([Rdiag Roffd]));
        
        Mlhat_fun = kernel.mollkernel_fourier(c, 0, 0, 0);
        Mlhat_fun1 = kernel.mollkernel_fourier(1, 0, 0, 0);
        %Mlhat_fun = kernel.diffkernel_fourier(c, 0, 0, 1);
        %Ftrunc(i) = max(abs(Mlhat_fun(f)));
        Ftrunc(i) = abs(kernel.fourier_scaling(c^2, 0))/c;
    end
    normalization = max(abs(kernel.direct([1 0 0], [0 0 0], f)));
    Rtrunc = Rtrunc / normalization;
end

function show_split(kernel, tol, style)
    level = 0;
    k = linspace(0, kernel.Kmax, 200);
    gamma_hat = kernel.fourier_scaling(k.^2, level);

    r = linspace(0.001, 1, 200);
    [Rdiag, Roffd] = kernel.real_decay(r, level);
    
    subplot(1,2,1)
    semilogy(k, abs(gamma_hat), style, displayname=kernel.name)
    hold on
    grid on
    ylim([-inf 1])
    xlabel('k')
    title('Fourier scaling')
            
    subplot(1,2,2)
    plot(r, abs(Roffd)+abs(Rdiag), style, displayname=kernel.name)
    hold on
    set(gca, yscale='log')
    grid on
    ylim([tol/1e2 1])
    xlabel('r')
    title('Real scaling')
end
