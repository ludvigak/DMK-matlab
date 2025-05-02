clear

set(groot, 'defaultTextInterpreter','latex');

N = 5000;
[filepath,~,~] = fileparts(mfilename('fullpath'));
datafile = sprintf("%s/data/parameter_data_N%d.mat", filepath, N);

ds = load(datafile, 'data_dict');
data_dict = ds.data_dict;

% Pick kernel
periodic = false;
if true
    name = kernels.stokeslet_pswf.name;
    shape = 'c';
else
    name = kernels.stokeslet_hasimoto.name;
    shape = 'sigma';
end
data = data_dict{name};
data =  data(data.periodic==periodic, :);
data.relerr_l2 = data.err_l2 ./ data.u_l2;



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
    semilogy(rescale(sub.(shape)), sub.relerr_l2, '-', LineWidth=2, displayname=sprintf('p=%d', p))
    hold on
end
grid on
xlabel(shape)
ylabel(error_label)
title('$p = 5:5:60$', interpreter='latex')
axis([0 inf 1e-15 1])
setup_fig();

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
    semilogy(sub.p, sub.relerr_l2, LineWidth=2)
    hold on
end
grid on
xlabel('$p$')
ylabel(error_label)
title('$c = (5:5:50) \pi/3$', interpreter='latex')
axis([0 60 1e-15 1])
setup_fig();


sfigure(3);
clf
x = rescale(best_p_data.(shape));
y = best_p_data.p;
plot(x, y, '.k')
%axis equal
fit_mask = best_p_data.relerr_l2 > 1e-13;
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
setup_fig();

sfigure(4);
clf
z = best_p_data.relerr_l2;
semilogy(x, z, '.')
P2 = polyfit(x(fit_mask), log(z(fit_mask)), 2)
hold on
semilogy(xx, exp( polyval(P2, xx) ), '-')
grid on
xlabel(xlab)
ylabel(error_label)
axis([0 inf 1e-16 1])
setup_fig();

sfigure(5);
clf
semilogy(best_p_data.p, best_p_data.relerr_l2, '.')
xlabel('$p$')
ylabel(error_label)
grid on
setup_fig();

write_fig(1, 'fig/stokeslet_pswf_c_errors')
write_fig(2, 'fig/stokeslet_pswf_p_errors')
write_fig(3, 'fig/stokeslet_pswf_c_vs_p')
write_fig(4, 'fig/stokeslet_pswf_c_err_fit')
write_fig(5, 'fig/stokeslet_pswf_p_err_fit')
