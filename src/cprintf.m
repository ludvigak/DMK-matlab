function cprintf(opt, varargin)
% Conditional printf
    if opt.verbose
        fprintf(varargin{:});
    end
end
