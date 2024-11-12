function u = dmk_apply(charges, dmk_state)
    opt = dmk_state.opt;
    % Unpack some opts
    tree = dmk_state.tree;
    p = dmk_state.opt.p;
    % Upward pass
    tic_up = tic();
    proxy_charges = init_proxy_charges(tree, charges, p, dmk_state.opt.kernel);
    t_up = toc(tic_up);
    cprintf(opt, '[dmk_apply] upward pass:        %.3f\n', t_up);
    % Downward pass
    tic_down = tic();
    local_expansions = collect_local_expansions(tree, ...
                                                proxy_charges, ...
                                                dmk_state.Twin, ...
                                                dmk_state.Tprox2pw, ...
                                                dmk_state.Tpw2poly, ...
                                                dmk_state.Tpwshift, ...
                                                dmk_state.Tp2c);
    t_down = toc(tic_down);
    cprintf(opt, '[dmk_apply] downward pass:      %.3f\n', t_down);
    % Local eval
    tic_eval = tic();
    ufar = evaluate_local_expansions(local_expansions, tree, p);
    t_eval = toc(tic_eval);
    cprintf(opt, '[dmk_apply] eval expansions:    %.3f\n', t_eval);
    % Near interactions
    tic_local = tic();
    ures = local_interactions(tree, charges, dmk_state.opt.kernel);
    t_local = toc(tic_local);
    cprintf(opt, '[dmk_apply] local interactions: %.3f\n', t_local);
    % Self interactions
    uself = self_interaction(tree, charges, dmk_state.opt.kernel);
    % Sum up
    u = ufar + ures + uself;
    cprintf(opt, '[dmk_apply]               SUM = %.3f\n', t_down+t_up+t_eval+t_local);        
end

