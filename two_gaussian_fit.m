function [ proportions,means,sigmas,mixture_model ] = ...
    two_gaussian_fit( velocities,V_high )
%TWO_GAUSSIAN_FIT Tries to fit a gaussian mixture model of two gaussians to
%   velocity data points
%
%   velocities - vector containing velocity data points
%   V_high - velocity up to which velocity data points are considered
%   proportions - proportional weight of the two Gaussians
%   means - means of the two Gaussians
%   variances - variances of the two Gaussians

velocities = real(velocities);
velocities = velocities(velocities<=V_high);

% Fit two Gaussian mixture model
fit_options = statset('Display','off','MaxIter',10000, ...
    'TolFun',1e-8);
mixture_fit = gmdistribution.fit( ...
    velocities.',2,'Options',fit_options);
mixture_model = @(vv) pdf(mixture_fit,vv);

% Two Gaussians mixture model parameters
proportions = mixture_fit.PComponents;
means = mixture_fit.mu;
sigmas = squeeze(mixture_fit.Sigma).';

% Sort parameters so that first component is that with lower mean
[means,sort_inds] = sort(means,'ascend');
proportions = proportions(sort_inds);
sigmas = sigmas(sort_inds);

end

