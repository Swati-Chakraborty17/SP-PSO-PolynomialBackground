%% -----------------------------
%% Objective Function
%% -----------------------------
function Jval = objective_function(params, x, V_obs, xlmt)

%   Jval = OBJECTIVE_FUNCTION(params, x, V_obs, xlmt) evaluates the
%   percentage RMS misfit between observed and modeled SP data and includes
%   an optional soft boundary penalty to prevent boundary-hugging artifacts 
%   in PSO to discourage solutions near parameter limits.
%
%   Inputs:
%       params : Model parameter vector
%       x      : Observation locations
%       V_obs  : Observed SP data
%       xlmt   : Parameter bounds [lower, upper]
%
%   Output:
%       Jval   : Objective function value 

    % Forward model response (thin-sheet source + polynomial background)
    V_calc = model_thin_sheet(x, params);

    % Percentage RMS misfit (normalized by observed data)
    misfit = 100 * sqrt(sum(((V_obs - V_calc) ./ V_obs).^2) / length(V_obs) );

    % Soft boundary penalty:
    % Introduced to discourage boundary-hugging solutions in PSO while
    % maintaining a smooth objective landscape.
    margin = 0.1;          % Fractional buffer from parameter bounds
    w = 0;                 % Penalty weight (set to zero if inactive)

    % Parameter ranges and effective interior bounds
    prange = xlmt(:,2) - xlmt(:,1);
    low  = xlmt(:,1) + margin * prange;
    high = xlmt(:,2) - margin * prange;

    % Compute penalty for parameters outside interior bounds
    penalty = 0;
    for i = 1:length(params)
        if params(i) < low(i)
            penalty = penalty + ((low(i) - params(i)) / prange(i))^2;
        elseif params(i) > high(i)
            penalty = penalty + ((params(i) - high(i)) / prange(i))^2;
        end
    end

    % Total objective value
    Jval = misfit + w * penalty;

end