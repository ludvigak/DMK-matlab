function dmk_opt = dmk_default_opts(args)
    arguments
        args.tolerance % tolerance
        args.kernel function_handle = @kernels.laplace_ewald
                       % which kernel+split to use
        args.p = -1    % polynomial order, default (-1) is to adaptively decide
        args.verbose (1,1) logical = false
    end
    % Unpack args
    p      = args.p;
    tol    = args.tolerance;
    kernel = args.kernel();
    % Setup mollification width
    sigma_0 = 1/sqrt(log(1/tol));
    if p==-1
        % Estimate suitable p
        p_list = 2:2:100;
        interp_err = estimate_interp_error(p_list, sigma_0, kernel);
        idx = find(interp_err < tol, 1);
        p = p_list(idx);
        p_err = interp_err(idx);
        if isempty(idx)
            error("Unable to find p satisfying tolerance")
        end
        cprintf(args, "[dmk_default_opts] selecting p=%d\n", p);
    else
        p_err = estimate_interp_error(p, sigma_0);
    end
    % Setup root level windowed kernel
    Kmax_win = ceil( 2*log(1/tol));
    nf_win = Kmax_win;
    hf_win = Kmax_win/nf_win;
    Ctrunc = sqrt(3) + 6*sigma_0;    
    % Setup planewave ops
    r0 = 1;
    D = 3*r0;
    h0 = 2*pi/D;
    K0 = 4/1 * log(1/tol);
    nf = ceil(K0/h0);
    % Store in struct
    dmk_opt = struct();
    dmk_opt.kernel = kernel;
    dmk_opt.p = p;
    dmk_opt.tol = tol;
    dmk_opt.sigma_0 = sigma_0;
    dmk_opt.nf_win = nf_win;
    dmk_opt.hf_win = hf_win;
    dmk_opt.Ctrunc = Ctrunc;
    dmk_opt.h0 = h0;
    dmk_opt.nf = nf;
    dmk_opt.verbose = args.verbose;
    % Report
    if dmk_opt.verbose
        opt_fields = fields(dmk_opt);
        for i=1:numel(opt_fields)
            name = opt_fields{i};
            if strcmp(name, "kernel")
                val = dmk_opt.kernel.name;
            else
                val  = num2str(dmk_opt.(name));
            end
            fprintf("[dmk_default_opts] %s = %s\n", name, val);
        end
    end
    % Print error estimates
    cprintf(dmk_opt, "[dmk_default_opts] estimated interp. rel.err=%.2e\n", p_err);
    trunc_err = norm(...
        kernel.diffkernel([1 0 0], [0 0 0], ones(1, kernel.dim_in), sigma_0), ...
        inf);
    cprintf(dmk_opt, "[dmk_default_opts] estimated kern. trunc.err=%.2e\n", trunc_err);
end
