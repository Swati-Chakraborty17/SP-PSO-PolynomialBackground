%This script performs opttimization of self-potential (SP) data using a
% dipping thin-sheet model superimposed on a quadratic regional background.
%
% The optimization is carried out using Particle Swarm Optimization (PSO)
% algorithm with multiple independent runs to assess solution robustness
% and stability.
%
% The background potential is modeled as:
%   V_bg(x) = c + b*x + i*x^2
%
% Model parameters:
%   k      : polarization strength (mV)
%   x0     : horizontal position of sheet center (m)
%   h      : depth to sheet center (m)
%   a      : half-length of the sheet (m)
%   alpha  : dip angle of the sheet (radians)
%   c, b, i: coefficients of quadratic background potential
%
% The PSO implementation is adapted and modified from standard PSO
% formulations (see comments in pso.m function for attribution).
%
% Author: Swati Chakraborty
% Department of Geology & Geophysics, IIT Kharagpur
% Date  : 31/03/2026
%% ------------------------------------------------------------------------

clc; clear; close all;
t_start = tic;     % Start timer

%% -----------------------------
%% Step 1: Read Observed SP Data
%% -----------------------------
data1 = readmatrix('Model_eta.dat');   % Load numeric matrix
x = data1(:,1);                        % Profile distance
V_obs =data1(:,2);                     % Observed values

%% -----------------------------
%% Step 2: PSO Setup
%% -----------------------------
nParams    = 8;      % Number of unknown model parameters
nParticles = 50;     % Swarm size
T          = 200;    % Number of PSO iterations
nRuns      = 10;     % Independent PSO runs for stability analysis

% Parameter search bounds for 'Model_eta.dat'

xlmt = [   370   , 570   ;    
           450   , 550   ;   
           15    , 120   ;    
           5     ,  80   ;    
    deg2rad(0)   , deg2rad(90); 
           0     ,  40   ;    
          -0.5   , 0     ;    
            0    , 1e-4 ];   

% Objective function
J = @(params) objective_function(params, x, V_obs, xlmt);

% Storage variables
X_all_good   = [];                      % All acceptable models
x_best_all   = zeros(nParams, nRuns);   % Best model per run
f_best_all   = zeros(1, nRuns);         % Best misfit per run
fitness_all  = zeros(nRuns, T+1);       % Convergence history

threshold    = 2;                       % threshold misfit
param_labels = {'k','x0', 'h', 'a','alpha (deg)', 'c', 'b' ,'i'};

%% -----------------------------
%% Step 3: Multiple PSO Runs
%% -----------------------------
for run = 1:nRuns
    
    [f_best, x_best, X_hist, F_hist, fg, g] = pso(J, nParams, xlmt, nParticles, T);
    x_best_all(:, run) = x_best;
    f_best_all(run)    = f_best;
    fitness_all(run,:) = fg;
    g_best_run = g;          % keep convergence history of parameters

    % Collect all models with acceptable misfit
    X_temp = reshape(X_hist, nParams, []);
    F_temp = reshape(F_hist, 1, []);
    good_idx = F_temp < threshold;
    X_all_good = [X_all_good, X_temp(:, good_idx)];

end

% Identify best run
[~, best_run_idx] = min(f_best_all);
fg_best_run = fitness_all(best_run_idx, :);

%% ----------------------------------------------
%% Step 4: Filter models based on h > a*sin(alpha)
%% ----------------------------------------------
% Enforces geometric constraint ensuring the sheet does not intersect the surface
h_all     = X_all_good(3,:);      % h
a_all     = X_all_good(4,:);      % a
alpha_all = X_all_good(5,:);      % alpha (in radians)

% Condition: h > a*sin(alpha)
valid_idx = h_all > (a_all .* sin(alpha_all));

% Keep only valid models
X_all_good = X_all_good(:, valid_idx);

%% ----------------------------------------------
%% Step 5: Print Results from All Runs
%% ----------------------------------------------
fprintf('\n--- Results from All Runs ---\n');
for r = 1:nRuns
    x_run = x_best_all(:, r);
    x_run(5) = rad2deg(x_run(5)); % convert alpha
    fprintf('\nRun %d:\n', r);
    for i = 1:nParams
        fprintf('%s = %.4f\n', param_labels{i}, x_run(i));
    end
    fprintf('RMS Misfit = %.4f%%\n', f_best_all(r));
end

% Best run summary
[~, best_run_idx] = min(f_best_all);
fprintf('\n--- Best Run = %d ---\n', best_run_idx);
x_best_deg = x_best_all(:, best_run_idx);
x_best_deg(5) = rad2deg(x_best_deg(5));
for i = 1:nParams
    fprintf('%s = %.4f\n', param_labels{i}, x_best_deg(i));
end
fprintf('RMS Misfit = %.4f%%\n', f_best_all(best_run_idx));

runtime = toc(t_start);   % Stop timer
fprintf('Execution time: %.2f seconds\n', runtime);

%% ------------------------------------------------------------
%% Step 6: Plot Parameter + Misfit Convergence Over Generations 
%% ------------------------------------------------------------
numIter = size(X_hist, 3);   % Total number of PSO generations stored in the swarm history
param_evolution = zeros(nParams, numIter);

for t = 1:numIter
    % param_evolution(:, t) = X_hist(:, best_run_idx, t);
    [~, best_idx] = min(F_hist(:, t));
    param_evolution(:, t) = X_hist(:, best_idx, t);
end

% Convert alpha to degrees for readability
param_evolution(5, :) = rad2deg(param_evolution(5, :));

% LaTeX-style parameter labels
param_labels = { ...
    'k (mV)', ...
    'x_{0} (m)', ...
    'h (m)', ...
    'a (m) ', ...
    '\alpha( \circ)',...
    'c (m)',...
    'b (mV / m)',...
    'i (mV/m^2)'
};

%% Plot :1 [Parameter+Misfit evolution]
figure('Color','w');

for i = 1:nParams + 1   % add one extra loop for misfit
    if i <= nParams
        % ----- Parameter convergence subplots -----
        ax = subplot(nParams+1, 3, (i-1)*3 + 1);
        plot(1:numIter, param_evolution(i,:), 'b-', 'LineWidth', 1.5);
        ylabel(param_labels{i}, 'Interpreter', 'tex', ...
               'FontSize', 12, 'FontWeight', 'bold');
    else
        % ----- Misfit reduction subplot -----
        ax = subplot(nParams+1, 3, (i-1)*3 + 1);
        plot(0:T, fg_best_run, 'b-', 'LineWidth', 1.5);
        ylabel('Misfit (%)', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Generations', 'FontSize', 12, 'FontWeight', 'bold');
        %title('Misfit Reduction', 'FontWeight', 'bold', 'FontSize', 12);
    end

    % ----- Common formatting for all subplots -----
    xlim([0 T]);
    ytickformat('%.1f');
    set(gca, 'FontWeight', 'bold', 'FontSize', 8, 'LineWidth', 1.2, 'Box', 'on');

    % Adjust subplot size (keeps all aligned)
    pos = get(ax, 'Position');
    pos(3) = 0.22;   % width
    pos(4) = 0.07;   % height
    set(ax, 'Position', pos);

    % Remove X tick labels for all except last (misfit)
    if i < nParams + 1
        set(gca, 'XTickLabel', []);
    end
end
ylim([0 100]);


%% -----------------------------------
%% Step 7: Statistical Analysis
%% Mean and Std of Good Models
%% -----------------------------------
N = size(X_all_good,2);
mean_model = mean(X_all_good,2);
std_dev    = std(X_all_good,0,2);

% Convert alpha to deg
mean_model_deg = mean_model;
mean_model_deg(5) = rad2deg(mean_model_deg(5));
std_dev_deg = std_dev;
std_dev_deg(5) = rad2deg(std_dev(5));

param_labels = {'k','x0', 'h', 'a','alpha (deg)','c', 'b','i'};

fprintf('\n--- Mean Model averaged from across %d Runs ---\n',nRuns);
for i = 1:nParams
    fprintf('%s = %.4f ± %.4f\n', param_labels{i}, mean_model_deg(i),std_dev_deg(i));
end

%% Plot: 2 [Adaptive Histogram distribution + Density Curve]

figure('Color','w');  % white background

% Parameter labels
param_labels = { ...
    'k (mV)', ...
    'x_{0} (m)', ...
    'h (m)', ...
    'a (m)', ...
    '\alpha( \circ)', ...
    'c (m)', ...
    'b (mV / m)', ...
    'i (mV/m^2)'
};

for i = 1:nParams
    subplot(3, 8, i);

    % Select data and statistics
    if i == 5
        data = rad2deg(X_all_good(i,:));
        mu = rad2deg(mean_model(i));
        sigma = rad2deg(std_dev(i));
    else
        data = X_all_good(i,:);
        mu = mean_model(i);
        sigma = std_dev(i);
    end

    % --- Adaptive bin width control ---
    n = numel(data);
    IQR_val = iqr(data);
    bin_width = 2 * IQR_val / (n^(1/3));  % Freedman–Diaconis base width

    % Compute number of bins and enforce limits 
    num_bins = round((max(data) - min(data)) / bin_width);
    num_bins = min(max(num_bins, 3),13);  % constrain for visibility

    edges = linspace(min(data), max(data), num_bins);

    % --- Histogram plotting ---
    h = histogram(data, edges, 'FaceColor', [0.4 0.7 0.9], ...
        'EdgeColor', 'k', 'LineWidth', 1.0, 'Normalization', 'count'); 
    hold on;

    % --- Density curve overlay (scaled to histogram) ---
    x_vals = linspace(min(data), max(data), 300);
    pdf_vals = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x_vals - mu)/sigma).^2);
    pdf_scaled = pdf_vals * max(h.Values) / max(pdf_vals);
    plot(x_vals, pdf_scaled, 'r-', 'LineWidth', 1.3);

    % --- Labels and formatting ---
    xlabel(param_labels{i}, 'Interpreter','tex','FontWeight','bold','FontSize',12);
    if i == 1 
        ylabel('Number of Models', 'FontWeight','bold', 'FontSize',12);
    end
    ax = gca;
    ax.XRuler.Exponent = 0;  % hide exponent from axis
    xt = ax.XTick;           % extract tick values
    
    ax.XLabel.Position(2) = ax.XLabel.Position(2) - 0.05*(ax.YLim(2)-ax.YLim(1)); % move down a bit

    ax = gca;
    ax.FontSize = 11;
    ax.LineWidth = 1.2;
    ax.TickDir = 'in';
    ax.FontWeight = 'bold';
    box off;
   
end

%% --------------------------------
%% Step 8: Correlation & Covariance
%% --------------------------------
cov_matrix = cov(X_all_good');
corr_matrix = corrcoef(X_all_good');

disp('Covariance Matrix:');
disp(cov_matrix);
disp('Correlation Matrix:');
disp(corr_matrix);

%% ------------------------------------------------------
%% Step:9 Mean Model within ±1 Standard Deviation
%% ------------------------------------------------------
% Compute mean and std of current X_all_good
mean_params = mean(X_all_good, 2);
std_params  = std(X_all_good, 0, 2);

% Logical index for models within ±1 std dev for all parameters
within_std_idx = all(abs(X_all_good - mean_params) <= std_params, 1);

X_within_std = X_all_good(:, within_std_idx);
fprintf('Number of models within ±1 std dev: %d\n', size(X_within_std,2));

% Compute mean model from this subset
mean_model_1std = mean(X_within_std, 2);
std_model_1std  = std(X_within_std, 0, 2);

% Convert alpha to degrees 
mean_model_1std_deg = mean_model_1std;
mean_model_1std_deg(5) = rad2deg(mean_model_1std_deg(5));

std_model_1std_deg = std_model_1std;
std_model_1std_deg(5) = rad2deg(std_model_1std_deg(5));

%% Display

fprintf('\n--- Mean Model (within 1 std dev) ---\n');
param_labels_deg = {'k','x0','h','a','alpha (deg)','c','b','i'};  
for i = 1:length(param_labels_deg)
    fprintf('%s = %.4f ± %.4f\n', ...
        param_labels_deg{i}, mean_model_1std_deg(i), std_model_1std_deg(i));
end

%% --------------------------------------
%% Step 10: Calculating model data 
%% --------------------------------------
% mean model obtained by averaging the best solutions from all independent PSO runs
V_best_10runs = model_thin_sheet(x, mean_model);
rms_misfit_best_10_runs = 100 * sqrt(sum(((V_obs - V_best_10runs)./V_obs).^2)/length(V_obs));
fprintf('\nRMS Misfit (Mean Model for models from 10 best runs) = %.4f%%\n', rms_misfit_best_10_runs);

% ensemble-derived mean model computed from solutions within ±1 standard deviation
V_mean = model_thin_sheet(x, mean_model_1std);
rms_misfit = 100 * sqrt(sum(((V_obs - V_mean)./V_obs).^2)/length(V_obs));
fprintf('\nRMS Misfit (Mean Model for models within 1 standard deviation) = %.4f%%\n', rms_misfit);

% Compute estimated background
for i =1: length(x)
    V_bg(i) = mean_model_1std_deg(6) + mean_model_1std_deg(7)*x(i) +mean_model_1std_deg(8)*x(i).^2 ;
end
figure('Color','w'); % white background

%% Plot: 3 [ Fitting between observed and calculated data]

plot(x, V_obs, 'o', 'MarkerSize', 5, 'LineWidth', 1.5, ...
     'MarkerFaceColor', [0 0 0.6], ...   % Deep blue fill
     'MarkerEdgeColor', [0 0 0.6]);hold on;      % Deep blue edge
%plot(x, V_mean_withoutbackground, 'k--', 'LineWidth', 1.7); hold on; %import from dat. file

% Plot mean model
plot(x, V_mean, 'r-', 'LineWidth', 1.7);hold on;
plot(x, V_bg, '-d', 'Color', [0 0.5 0], ...
     'LineWidth', 1.7, 'MarkerSize', 2, ...
     'MarkerFaceColor', [0 0.5 0]);

% Legend
legend({'Observed data','Model data with background','Background'}, 'Location','northeast', ...
       'FontSize',12,'FontWeight','bold', 'Interpreter','tex', 'Box','off');

% Labels
xlabel('Distance (m)', 'FontSize',14,'FontWeight','bold');
ylabel('S.P. Anomaly (mV)', 'FontSize',14,'FontWeight','bold');

% Grid and axis formatting
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.TickDir = 'in';
ax.Box = 'on';
% Make tick labels bold
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

% Tight layout
xlim([min(x) max(x)]);

%% Plot: 4 [thin-sheet geometry] 

figure('Color','w'); hold on;

% Extract mean parameters
k     = mean_model_1std(1);
x0    = mean_model_1std(2);   % center (horizontal)
h     = mean_model_1std(3);   % depth to center
a     = mean_model_1std(4);   % half-length
alpha = mean_model_1std(5);   % dip angle (rad)

% Direction vector (alpha from +x axis)
dx = cos(alpha);
dz = sin(alpha);

% Endpoints symmetric about center
x1 = x0 + a*dx;
z1 = h - a*dz;

x2 = x0 - a*dx;
z2 = h + a*dz;

% Plot sheet
plot([x1 x2],[z1 z2],'-','Color',[0 0.6 0],'LineWidth',7); 

% Extend vertical line from sheet center to surface
plot([x0 x0],[0 h],'k--','LineWidth',1.2);

% Horizontal line from left axis (xl(1)) to sheet center
plot([-160 x0],[h h],'k--','LineWidth',1.2);

% Grid and axis formatting
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.TickDir = 'in';
ax.Box = 'on';

% Make tick labels bold
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

% Mark sheet center
plot(x0,h,'go','MarkerFaceColor','k','MarkerSize',6);

% Axes formatting
set(gca,'YDir','reverse'); % depth downwards
xlabel('Distance (m)','FontSize',13,'FontWeight','bold');
ylabel('Depth (m)','FontSize',13,'FontWeight','bold');
axis equal; box on;
ylim([0 max([z1 z2])*1.6]);
xlim([min(x)-20 max(x)+20]);

% Get current axis limits
xl = xlim;
yl = ylim;
% Place parameter text in top-right corner
text(xl(2)-0.05*range(xl), yl(1)+0.5*range(yl), ...
    {sprintf('k = %.2f mV',k),...
     sprintf('x_0 = %.2f m',x0), ...
     sprintf('h = %.2f m',h), ...
     sprintf('a = %.2f m',a), ...
     sprintf('\\alpha = %.2f°',rad2deg(alpha))}, ...
     'FontSize',11,'FontWeight','bold','BackgroundColor','none','EdgeColor','none', ...
     'HorizontalAlignment','right');

title('Estimated Thin Sheet Geometry','FontSize',14 ,'FontWeight','bold');

%% Plot: 5 [Background Potential, V_bg]

figure('Color','w'); % white background
plot(x, V_bg, '-d', 'Color', [0 0.5 0], ...
     'LineWidth', 1.7, 'MarkerSize', 2, ...
     'MarkerFaceColor', [0 0.5 0]);
% Labels
xlabel('Distance (m)', 'FontSize',14,'FontWeight','bold');
ylabel('S.P. Anomaly (mV)', 'FontSize',14,'FontWeight','bold');
title('Estimated Background Variation','FontSize',14 ,'FontWeight','bold');

%% Plot: 6 [Conditional Cost Maps: Source vs Background]

param_labels = { ...
    'k (mV)', ...
    'x_{0} (m)', ...
    'h (m)', ...
    'a (m) ', ...
    '\alpha( \circ)',...
    'c (m)',...
    'b (mV / m)',...
    'i (mV/m^2)'
};

source_params_idx = [1, 2, 3, 4, 5]; % k, x0, h, a, alpha
bg_params_idx     = [6, 7, 8];       % c, b, i

nGrid = 50;
figure('Color','w');

for r = 1:length(bg_params_idx)
    for c_idx = 1:length(source_params_idx)

        src_idx = source_params_idx(c_idx);
        bg_idx  = bg_params_idx(r);

        % Grid (convert alpha to degrees ONLY for plotting)
        if src_idx == 5
            src_vals = linspace(rad2deg(xlmt(src_idx,1)), ...
                                rad2deg(xlmt(src_idx,2)), nGrid);
        else
            src_vals = linspace(xlmt(src_idx,1), xlmt(src_idx,2), nGrid);
        end
        bg_vals  = linspace(xlmt(bg_idx,1), xlmt(bg_idx,2), nGrid);

        Z = nan(nGrid);

        for i = 1:nGrid
            for j = 1:nGrid

                params = mean_model;

                % Assign values (convert alpha back to radians)
                if src_idx == 5
                    params(src_idx) = deg2rad(src_vals(i));
                else
                    params(src_idx) = src_vals(i);
                end
                params(bg_idx) = bg_vals(j);

                fixed_idx = [src_idx, bg_idx,2,3,4];
                free_idx  = setdiff(1:length(params), fixed_idx);

                if ~isempty(free_idx)
                    fun = @(p) objective_function( ...
                        update_params(params, free_idx, p), ...
                        x, V_obs, xlmt);
                    opt_vals = fminsearch(fun, params(free_idx), ...
                        optimset('Display','off'));
                    Z(j,i) = fun(opt_vals);
                else
                    Z(j,i) = objective_function(params, x, V_obs, xlmt);
                end
            end
        end

        % Smooth contours (visual only)
        [Xs, Ys] = meshgrid(src_vals, bg_vals);
        Zs = interp2(Xs, Ys, Z, Xs, Ys, 'cubic');

        subplot(length(bg_params_idx), length(source_params_idx), ...
            (r-1)*length(source_params_idx) + c_idx);

        % Filled contours
        contourf(src_vals, bg_vals, Zs, 40, 'LineColor','none');
        hold on;

        levels = linspace(min(Zs(:)), max(Zs(:)), 40);
        targetLevel = levels(3);
        
        contour(src_vals, bg_vals, Zs, ...
                [targetLevel targetLevel], ...
                'k--', 'LineWidth', 1.5);


        set(gca,'YDir','normal');
        xlabel(['$', param_labels{src_idx}, '$'], ...
            'Interpreter','latex');
        ylabel(['$', param_labels{bg_idx}, '$'], ...
            'Interpreter','latex');

        % Mark mean solution
        if src_idx == 5
            plot(rad2deg(mean_model(src_idx)), mean_model(bg_idx), ...
                'wo','MarkerSize',4,'LineWidth',1.5);
        else
            plot(mean_model(src_idx), mean_model(bg_idx), ...
                'wo','MarkerSize',4,'LineWidth',1.5);
        end

        axis tight;
        hold off;
    end
end

sgtitle('Conditional Cost Maps: Source vs Background Parameters', ...
    'FontWeight','bold','FontSize',14);

% Helper function to update parameters
function new_params = update_params(params, idx, values)
    new_params = params;
    new_params(idx) = values;
end

%% ------------------------------------------------
%% Step 11: Save the results as an Output .mat file
%% ------------------------------------------------
save('alpha_results.mat', ...
    'V_obs','xlmt','fg_best_run','X_all_good','param_evolution', ...
    'mean_model','std_dev','X_within_std','mean_model_1std','std_model_1std', ...
    'V_best_10runs','V_mean');