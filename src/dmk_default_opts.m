function dmk_opt = dmk_default_opts(args)
    arguments
        args.tolerance % tolerance
        args.kernel function_handle = @kernels.laplace_ewald
                       % which kernel+split to use
        args.p = -1    % polynomial order, default (-1) is to adaptively decide
        args.periodic (1,1) logical = false;
        args.verbose (1,1) logical = false
    end
    % Unpack args
    p      = args.p;
    tol    = args.tolerance;
    kernel = args.kernel(tolerance=tol);
    if p==-1
        % Estimate suitable p
        p_list = 2:2:100;
        interp_err = estimate_interp_error(p_list, kernel);
        idx = find(interp_err < tol, 1);
        p = p_list(idx);
        p_err = interp_err(idx);
        if isempty(idx)
            error("Unable to find p satisfying tolerance=%g, lowest interp error was %g", ...
                  tol, min(interp_err))
        end
        cprintf(args, "[dmk_default_opts] selecting p=%d\n", p);
    else
        p_err = estimate_interp_error(p, kernel);
    end
    % TODO: Move default parameters into kernel
    % Setup root level windowed kernel
    Kmax_win = ceil(kernel.Kmax);
    nf_win = Kmax_win;
    hf_win = Kmax_win/nf_win;
    Ctrunc = sqrt(3) + 1;
    % Root level periodic
    nf_per = ceil(kernel.Kmax / (2*pi));
    % Setup planewave ops
    r0 = 1;
    D = 3*r0;
    h0 = 2*pi/D;
    K0 = 2*kernel.Kmax;
    nf = ceil(K0/h0);
    % Store in struct
    dmk_opt = struct();
    dmk_opt.kernel = kernel;
    dmk_opt.p = p;
    dmk_opt.tol = tol;
    dmk_opt.nf_win = nf_win;
    dmk_opt.hf_win = hf_win;
    dmk_opt.Ctrunc = Ctrunc;
    dmk_opt.nf_per = nf_per;
    dmk_opt.h0 = h0;
    dmk_opt.nf = nf;
    dmk_opt.periodic = args.periodic;
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
        kernel.diffkernel([0.9999999 0 0], [0 0 0], ones(1, kernel.dim_in), 0), ...
        inf);
    cprintf(dmk_opt, "[dmk_default_opts] estimated kern. trunc.err=%.2e\n", trunc_err);
end
