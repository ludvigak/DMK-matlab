function dmk_state = dmk_init(points, max_level, dmk_opt)
    dmk_state = struct(opt=dmk_opt);
    % Unpack opts
    p = dmk_opt.p;
    sigma_0 = dmk_opt.sigma_0;
    h0 = dmk_opt.h0;
    nf = dmk_opt.nf;
    % Build tree
    atic = tic();
    dmk_state.tree = octree(points, max_level);
    t_tree = toc(atic);
    % Precompute operators
    atic = tic();
    dmk_state.Twin = operator_windowed(p, dmk_opt.hf_win, dmk_opt.nf_win, ...
                                       dmk_opt.Ctrunc, sigma_0, dmk_opt.kernel);
    dmk_state.Tprox2pw = operator_proxy2planewave(p, h0, nf, max_level);
    dmk_state.Tpw2poly = operator_planewave2local(p, h0, nf, max_level, sigma_0, ...
                                                 dmk_opt.kernel);
    dmk_state.Tpwshift = operator_planewave_shift(h0, nf);
    dmk_state.Tp2c = operator_parent2child(p);
    t_ops = toc(atic);
    % Report
    cprintf(dmk_opt, '[dmk_init] build tree:           %.3f s\n', t_tree);
    cprintf(dmk_opt, '[dmk_init] precompute operators: %.3f s\n', t_ops);
end

