xi = 8; % Ewald parameter-- affects convergence of the real & fourier sums
%% Inpur parameters
%-------------------------------------------
icase = 3;
if icase == 1
    N = 1;  % Number of particles
    eval_idx = 1;    % index of points at which Stokeslet is evaluated
    x = [0 0 0;];   % Stokeslet Locations
    f0 = 1;  % Source strength
 
    f = [f0 1e-16 1e-16];   % not neutrally charged
elseif icase == 2
    N = 2;  % Number of particles
    eval_idx = 1:2;    % index of points at which Stokeslet is evaluated
    x = [0 0 0; .5 0 0];   % Stokeslet Locations
    f0 = 1;  % Source strength
    f = [f0 1e-16 1e-16; -f0 1e-16 1e-16];  % neutrally charged
elseif icase == 3
    N = 3;  % Number of particles
    eval_idx = 1:3;    % index of points at which Stokeslet is evaluated
    x = [0 0 0; .5 0 0; 0.8 0 0];   % Stokeslet Locations
    f0 = 1;  % Source strength
    f = [f0 1e-16 1e-16; -f0 1e-16 1e-16; 0 0 0];   
    % neutrally charged. 
    % added a zero strength stokeslet at the same position as the first stokeslet
end

box = [1 1 1];
SE_opt.M = 128;
SE_opt.box = box;
SE_opt.P = 25;
SE_opt.sampling_factor=5;

[u,ufft,ureal] = stokeslet_total2P(eval_idx,x,f,xi,SE_opt);
u
ufft
ureal