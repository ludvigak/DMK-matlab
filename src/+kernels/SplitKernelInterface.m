classdef SplitKernelInterface < handle
% Abstract class for interfacing to split kernels
    properties (Abstract, Constant)
        name
        % Name of the split kernel
        dim_in
        % Dimension of kernel input data, eg. 1 for Laplace, 3 for Stokeslet
        dim_out
        % Dimension of kernel output data, eg. 1 for Laplace, 3 for Stokeslet
    end
    methods (Abstract)
        direct(targets, sources, data)
        % Direct evaluation of the kernel
        % targets: (Ntrg, 3)
        % sources: (Nsrc, 3)
        % data: (Nsrc, dim_in)
        % returns: (Ntrg, dim_out)

        reskernel
        
        diffkernel
        
        diffkernel_fourier
        
        winkernel_fourier

        self_interaction
    end
end
