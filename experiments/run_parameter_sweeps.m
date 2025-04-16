% Run parameter sweeps to generate data for plotting
clear

% ==========
% SETUP

% Test setup
tree_max_level = 1;
N = 1000;

datafile = sprintf("parameter_data_N%d.mat", N);
replace_data = true; % Overwrite existing parameters tests

p_list = 2:60;

periodic_status = [true false];

nf_list = 2:50;

% PSWF
c_list = nf_list * pi/3;
pswf_kernels = {@kernels.rotlet_pswf, @kernels.stokeslet_pswf, @kernels.stresslet_pswf};

% Gaussian
sigma_list = sqrt(6 ./ (pi*nf_list));
gauss_kernels = {@kernels.rotlet_ewald, @kernels.stokeslet_hasimoto, @kernels.stresslet_hasimoto};

% ===========

% Load or create datafile
if isfile(datafile)
    ds = load(datafile, 'data_dict');
    data_dict = ds.data_dict;
else
    data_dict = dictionary();
    data_dict{"foo"} = table();
    data_dict.remove("foo");
end

tic_start = tic();

% Start running all kernels
kernel_list = [pswf_kernels, gauss_kernels];
for kernel_idx = 1:numel(kernel_list)
    kernel = kernel_list{kernel_idx};
    kernel_inst = kernel(tolerance=1e-10); % Dummy instantiation of kernel (to get name)
    kernel_name = kernel_inst.name;
    pswf_kernel = contains(kernel_name, 'PSWF');
    if pswf_kernel
        parameter_list = c_list;
    else
        parameter_list = sigma_list;
    end
    % Setup system
    rng(1);
    points    = rand(N, 3)-1/2;
    charges   = rand(N, kernel_inst.dim_in)-1/2;
    if kernel_inst.dim_in==6
        % Stresslet
        q = rand(N, 3)-1/2;
        n = rand(N, 3)-1/2;
        n = n ./ sqrt( sum(n.^2, 2) );
        charges = [q n];
    end
    % Always try to make system charge neutral (not just when periodic)
    charges = charges - sum(charges, 1)/N;
    charges(end, :) = charges(end, :)-sum(charges, 1);
    assert(all(abs(sum(charges)) < 1e-14))

    if data_dict.isKey(kernel_name)
        kernel_data = data_dict{kernel_name};
    else
        kernel_data = table([], [], [], [], [], [], [], [], ...
                            VariableNames={'periodic', 'p', 'c', 'sigma', 'err_l2', 'err_inf', 'u_l2', 'u_inf'});
    end
    for periodic=periodic_status
        if periodic
            u_ref = ewald_sum(points, charges, kernel_inst);
        else
            u_ref = kernel_inst.direct(points, points, charges);
        end
        u_l2 = norm(u_ref(:), 2);
        u_inf = norm(u_ref(:), inf);
        for p=p_list
            for param=parameter_list
                fprintf('[%s] param=%05.2f, p=%02d, periodic=%d. ', kernel_name, param, p, periodic*1);
                atic = tic();
                if pswf_kernel
                    args = {'c_pswf', param};
                    c = param;
                    sigma = 0;
                else
                    args = {'sigma_0', param};
                    c = 0;
                    sigma = param;
                end
                % Check if run already in table, then replace
                data_idx = find(kernel_data.periodic==periodic & ...
                                kernel_data.p==p & ...
                                kernel_data.c==c & ...
                                kernel_data.sigma==sigma, ...
                               1);
                if isempty(data_idx)
                    % Insert add end
                    data_idx = size(kernel_data, 1) + 1;
                else
                    fprintf('Run exists, ');
                    if replace_data
                        % Replace at idx
                        fprintf('replacing. ')
                    else
                        fprintf('skipping. \n');
                        continue;
                    end
                end
                dummy_tol = 1e-16;
                dmk_opt   = dmk_default_opts(verbose=false, ...
                                             kernel=kernel, ...
                                             tolerance=dummy_tol, ...
                                             p=p, ...
                                             kernel_args=args, ...
                                             periodic=periodic);
                dmk_state = dmk_init(points, tree_max_level, dmk_opt);
                u_dmk     = dmk_apply(charges, dmk_state);
                err_l2 = norm(u_ref(:) - u_dmk(:));
                err_inf = norm(u_ref(:) - u_dmk(:), inf);
                % Insert data
                kernel_data(data_idx, :) = ...
                    table(periodic, p, c, sigma, err_l2, err_inf, u_l2, u_inf);
                toc(atic)
            end
            data_dict{kernel_name} = kernel_data;
            save(datafile, 'data_dict');
        end
    end
end
toc(tic_start);
