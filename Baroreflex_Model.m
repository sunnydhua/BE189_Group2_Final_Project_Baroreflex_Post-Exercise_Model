clear; clc; close all;
%% Baroreflex Model
% Model tested over various hydration and training values

% ==========================
% PARAMETERS
% ==========================

% Hydration level
% 1.0 = hydrated
% 0.98 = moderate dehydration
% 0.95 = severe dehydration
alpha_hyd_vals = [1 0.98 0.95];

% Athletic training sweep
alpha_train_vals = 0.5:0.1:1.5;

% ==========================
% STORAGE ARRAYS
% ==========================
nTrain = length(alpha_train_vals);

max_error_all   = zeros(nTrain,1);
latency_all     = zeros(nTrain,1);

latency_plot           = zeros(nTrain,3);
BRE_plot               = zeros(nTrain,3);

% ======================================================
% SWEEPS
% ======================================================

for hyd_loop = 1:3
    
    alpha_hyd = alpha_hyd_vals(hyd_loop);

    if hyd_loop == 1
    currentColor = [1 0 0];
    elseif hyd_loop == 2
        currentColor = [0.9 0 0.9];
    else
        currentColor = [0 0 1];
    end

    for train_loop = 1:nTrain
    
        alpha_train = alpha_train_vals(train_loop);

        [t, BRE, MAP_error, P_star, MAP] = ...
            Baroreceptor_Calculations(alpha_hyd, alpha_train);
        
        % Metrics for this alpha_train

        [max_error_all(train_loop),idx] = max(MAP_error((t>0.1)));
        
        latency_all(train_loop) = t(idx);

    
        % PLOTS IN LOOP

        figure(1);
        sgtitle('BRE With Increasing Alpha_{train}')
        subplot(3,1,hyd_loop)
        latency_plot(train_loop) = plot(t, BRE, '-', 'LineWidth', 2, Color=currentColor .* [train_loop * 0.09 0 train_loop * 0.09]); hold on
        ylabel('BRE = -dHR/dMAP (bpm/mmHg)')
        ylim([-90 90])
        title(['Alpha_{hydration} = ', num2str(alpha_hyd)])
        grid on
        xlabel(subplot(3,1,3),'Time post-exercise (min)')

    
        figure(2)
        sgtitle('Max Distance From P* With Increasing Alpha_{train}')
        subplot(3,1,hyd_loop)
        BRE_plot(train_loop)  = plot(t, MAP_error, '-', 'LineWidth', 2, Color=currentColor .* [train_loop * 0.09 0 train_loop * 0.09]); hold on
        ylabel('MAP error = P^* - MAP (mmHg)')
        ylim([0 20])
        title(['Alpha_{hydration} = ', num2str(alpha_hyd)])
        grid on
        xlabel(subplot(3,1,3),'Time post-exercise (min)')
    
        % Find time that pressure changes direction
        max_error = max(MAP_error);
        latency = t(MAP_error == max_error);
    
    end % End train loop
    
    
    %% GENERAL PLOTS
   
    % Latency Plot
    figure(3)
    hold on
    plot(alpha_train_vals,latency_all, '*-', LineWidth=3, Color = currentColor)
    xlabel("Alpha_{ train}")
    xlim([0.45 1.55])
    ylabel("Latency (minutes)")
    ylim([0 20])
    title("Latency Time vs Alpha_{train}")
    legend(latency_plot([1 6 11]), "Alpha_{train} = 0.5", ...
       "Alpha_{train} = 1.0", ...
       "Alpha_{train} = 1.5", ...
       'FontSize', 10, 'FontWeight', 'Bold')
    grid on
    
    % Error Plot
    figure(4)
    hold on
    plot(alpha_train_vals,max_error_all, '*-', LineWidth=3, Color = currentColor)
    xlabel("Alpha_{ train}")
    xlim([0.45 1.55])
    ylabel("Max Distance From P* (mmHg)")
    ylim([0 20])
    title("Distance from P* vs Alpha_{train}")
    legend(BRE_plot([1 6 11]), "Alpha_{train} = 0.5", ...
       "Alpha_{train} = 1.0", ...
       "Alpha_{train} = 1.5", ...
       'FontSize', 8, 'FontWeight', 'Bold')
    grid on

end % end hydration loop

figure(3)
legend("Alpha_{hydration} = 1", ...
       "Alpha_{hydration} = 0.98", ...
       "Alpha_{hydration} = 0.95", ...
       'FontSize', 12, 'FontWeight', 'Bold')

figure(4)
legend("Alpha_{hydration} = 1", ...
       "Alpha_{hydration} = 0.98", ...
       "Alpha_{hydration} = 0.95", ...
       'FontSize', 12, 'FontWeight', 'Bold')