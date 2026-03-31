%% ----------------------------------------------------------------
%% Forward Model: 2D Inclined Thin Sheet with Quadratic Background
%% ----------------------------------------------------------------
function V = model_thin_sheet(x, params)

    %% --------------------------
    % Source parameters (SP anomaly)
    % k     : electric dipole moment (mV)
    % xa    : horizontal position of sheet center (m)
    % h     : depth to sheet center (m)
    % a     : half-length of the thin sheet (m)
    % alpha : dip angle of the sheet (radians)
    %% --------------------------
    k     = params(1);
    xa    = params(2);
    h     = params(3);
    a     = params(4);
    alpha = params(5);

    %% --------------------------
    % Background parameters
    % c : constant (offset)
    % b : linear regional trend
    % i : quadratic curvature term
    %
    % The quadratic polynomial (c + b·x + i·x²) represents a smooth regional
    % background commonly arising from large-scale geological or topographic
    % effects that vary slowly relative to the localized SP anomaly.
    %% --------------------------
    c     = params(6);
    b     = params(7);
    i     = params(8);

    %% --------------------------
    % Forward response computation
    %% --------------------------
    V = zeros(size(x));
    for j = 1:length(x)
        xj = x(j);

        % Thin-sheet SP anomaly (logarithmic formulation)+ quadratic background
        
        V(j) = k * log(((xj - xa - a*cos(alpha))^2 + (h - a*sin(alpha))^2) / ...
            ((xj - xa + a*cos(alpha))^2 + (h + a*sin(alpha))^2) )+ c + b * xj + i * xj.^2;
    end
end