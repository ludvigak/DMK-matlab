% Plot comparison of PSWF/Gauss splits for Stresslet

clear

c = 32; % 1e-10
sigma = 0.192;

Tp = kernels.stresslet_pswf(c_pswf=c);
Tg = kernels.stresslet_hasimoto(sigma_0=sigma);

% Real space end values (should be ~same)
%TpEND = Tp.Roffd_cheb(1)
%[~, TgEND] = Tg.real_decay(1, 0)

linewidth=1;
markersize=8;
% Real space plots
r = linspace(0, 1);
for f=1:2
    sfigure(f);
    clf
    [Tgdiag, Tgoffd] = Tg.real_decay(r, 0);
    plot(r, Tp.Roffd_cheb(r), linewidth=linewidth)
    hold on
    plot(r, Tgoffd, linewidth=linewidth)
    if f==1
        set(gca, yscale='linear')
        ylim([0 1.1])
    else
        set(gca, yscale='log')
        ylim([1e-11 10])
        a = gca;
        a.YAxis.TickValues = 10.^-(10:-2:0);
        a.YAxis.MinorTickValues = [];
    end
    grid on
    legend('PSWF','Gaussian', interpreter='latex')
    ylabel('$|T_{\textrm{offd}}(r)|$', interpreter='Latex')
    xlabel('$r$', interpreter='Latex')
end

% Fourier space plot
sfigure(3);clf;
k = linspace(1, Tp.Kmax);
pfplot = semilogy(k, abs(Tp.fourier_scaling(k.^2, 0)), linewidth=linewidth);
hold on
k = linspace(1, Tg.Kmax);
gfplot = semilogy(k, abs(Tg.fourier_scaling(k.^2, 0)), linewidth=linewidth);
grid on
plot(c, abs(Tp.fourier_scaling(c^2, 0)), '.', color=pfplot.Color, markersize=markersize)
plot(Tg.Kmax, abs(Tg.fourier_scaling(Tg.Kmax^2, 0)), '.', color=gfplot.Color, markersize=markersize)
ylabel('$|\widehat\gamma(k)|$', interpreter='Latex')
xlabel('$k$', interpreter='Latex')
legend('PSWF','Gaussian', interpreter='latex')
ylim([1e-11 10])
a = gca;
a.YAxis.TickValues = 10.^-(10:-2:0);
a.YAxis.MinorTickValues = [];

Kmaxfactor = c/Tg.Kmax

[filepath,~,~] = fileparts(mfilename('fullpath'));
fsize = [7.3, 5.8] * 0.9;
for f=1:3
    sfigure(f);
    setup_fig(width_cm=fsize(1), height_cm=fsize(2));
    write_fig(f, [filepath '/fig/split_compare_stresslet_' num2str(f)])
end
