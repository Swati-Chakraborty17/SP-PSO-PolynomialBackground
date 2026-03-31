%% ---------------------------------------------------------
%% PSO Function
%% ---------------------------------------------------------
% This Particle Swarm Optimization (PSO) implementation is
% based on standard PSO formulations described in:
% Elkmay (2026). pso - Particle Swarm Optimization (https://github.com/ElkmanY/pso), GitHub.

% The original structure has been modified for the present
% geophysical inverse problem by incorporating:
%   - Problem-specific parameter bounds
%   - Iteration-dependent linearly decreasing inertia weight
%   - Boundary reflection handling
%   - Storage of particle and fitness histories
%
% Author modifications are specific to SP anomaly optimization and background-aware modeling.
%% ---------------------------------------------------------

function [f_best, x_best, X_hist, F_hist, fg, g] = pso(J, d, xlmt, n, T)

    %% -----------------------------
    % Parameter bounds and velocity limits
    %% -----------------------------
    xb = xlmt(:,1);                % Lower bounds of parameters
    xp = xlmt(:,2);                % Upper bounds of parameters

    vb = -0.5*(xp-xb);             % Minimum velocity (relative to search range)
    vp =  0.5*(xp-xb);             % Maximum velocity

    %% -----------------------------
    % PSO control parameters
    %% -----------------------------
    c1 = 1.7;                      % Cognitive (personal learning) coefficient
    c2 = 1.8;                      % Social (global learning) coefficient
    w0 = 0.9;                      % Initial inertia weight
    w1 = 0.3;                      % Final inertia weight (linearly decreased)

    %% -----------------------------
    % Global best storage
    %% -----------------------------
    fg = zeros(1,T+1);             % Best fitness value per generation
    g  = zeros(d,T+1);             % Best solution vector per generation

    %% -----------------------------
    % Initialize particles
    %% -----------------------------
    X = rand(d,n).*(xp-xb) + xb;   % Particle positions uniformly within bounds
    V = rand(d,n).*(vp-vb) + vb;   % Initial particle velocities
    P = X;                         % Personal best positions

    %% -----------------------------
    % Evaluate initial fitness
    %% -----------------------------
    Fp = arrayfun(@(i) J(X(:,i)), 1:n);   % Fitness of all particles
    [fg(1), ig] = min(Fp);                % Global best fitness
    g(:,1) = X(:,ig);                     % Global best position

    %% -----------------------------
    % History storage (for convergence analysis)
    %% -----------------------------
    X_hist = zeros(d,n,T);          % Particle positions per iteration
    F_hist = zeros(n,T);            % Particle fitness per iteration

    %% -----------------------------
    % Main PSO loop
    %% -----------------------------
    for t = 1:T

        % Linearly decreasing inertia weight
        w = w0 - (w0 - w1) * (t / T);

        % Velocity update (inertia + cognitive + social terms)
        V = w*V ...
            + c1*(P - X).*rand(d,n) ...
            + c2*(g(:,t) - X).*rand(d,n);

        % Position update
        X = X + V;

        %% -----------------------------
        % Boundary handling (reflection)
        % Prevents particles from stagnating at bounds
        %% -----------------------------
        for j = 1:d
            X(j, X(j,:) < xb(j)) = xb(j) + 0.1*rand*(xp(j)-xb(j));
            X(j, X(j,:) > xp(j)) = xp(j) - 0.1*rand*(xp(j)-xb(j));
        end

        %% -----------------------------
        % Fitness evaluation
        %% -----------------------------
        F_ = arrayfun(@(i) J(X(:,i)), 1:n);

        % Update personal bests
        improved = F_ < Fp;
        P(:, improved) = X(:, improved);
        Fp(improved)   = F_(improved);

        % Update global best
        [fg(t+1), ig] = min(Fp);
        g(:, t+1) = P(:, ig);

        % Store history
        X_hist(:,:,t) = X;
        F_hist(:,t)   = F_;
    end

    %% -----------------------------
    % Final outputs
    %% -----------------------------
    x_best = g(:,end);              % Best solution found
    f_best = fg(end);               % Best fitness value
end