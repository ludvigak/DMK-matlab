% Development script for checking radial fourier kernels improved using Taylor expansions

clear all
clf
clc

xlog = logspace(-14, 0, 1000);
x = linspace(0, 1, 1000);

[f, df, d2f, d3f] = radial_fourier_kernels();
F ={f, df, d2f, d3f};

tic
[gf, gdf, gd2f, gd3f] = radial_fourier_kernels_gen();
toc
G = {gf, gdf, gd2f, gd3f};

k = [30]; % Will be O([1,50])

for D = 1:3
    sfigure(D);
    clf
    f = F{D+1};
    g = G{D+1};
    y0 = f(k, xlog);
    gdf(k, xlog);
    [y, err] = g(k, xlog);

    subplot(2,2,1)
    loglog(xlog, abs(y0), '.-', displayname='direct');
    hold on
    loglog(xlog, abs(y), '-', displayname='taylor');
    grid on
    title(sprintf('|d%df(k, r)|,   k=%d', D, k))
    legend()

    subplot(2,2,3)
    scale = max(abs(y),[],2);
    loglog(xlog, err ./ scale, '-', displayname='taylor accuracy');
    hold on
    loglog(xlog, eps + abs(y0-y) ./ scale, displayname='|taylor-direct|')
    loglog(xlog, eps*(1./xlog).^D * 10^(D-1) ./ k.^D, '-', displayname='error estimate')

    cutoff = (10^(D-1) ./ k.^D).^(1/D)
    loglog(cutoff, eps, '.', markersize=10, displayname='cutoff')

    grid on
    legend(location='SouthWest')
    ylim([1e-20, 1])

    subplot(2,2,2)
    plot(x, d3f(k, x))
    hold on
    plot(x, gd3f(k, x),'.')
    grid on
    title(sprintf('d%df(k, r),   k=%d', D, k))
end


