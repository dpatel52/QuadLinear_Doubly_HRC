function plot_moment_vs_curvature(Curvature, Envelope, beta_z1, beta_z2, beta_z3, beta_z4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Moment vs Curvature Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine the lengths of each zone
    len_beta_z1 = length(beta_z1);
    len_beta_z2 = length(beta_z2);
    len_beta_z3 = length(beta_z3);
    len_beta_z4 = length(beta_z4);

    % Calculate the indices at which each zone ends
    end_z1 = len_beta_z1;
    end_z2 = end_z1 + len_beta_z2;
    end_z3 = end_z2 + len_beta_z3;
    end_z4 = end_z3 + len_beta_z4;

    % Combine the indices into a vector
    beta_all2 = [end_z1; end_z2; end_z3; end_z4];

    % Define colors for each zone
    colors = {'k', 'm', 'b', 'r'};

    % Define labels for the legend
    zone_labels = {'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4'};

    % Plot each zone with a different color
    figure;
    hold on;
    for i = 1:length(beta_all2)
        if i == 1
            start_idx = 1;
        else
            start_idx = beta_all2(i-1) + 1;
        end
        end_idx = beta_all2(i);
        plot(Curvature(start_idx:end_idx), Envelope(start_idx:end_idx, 2), '-', 'LineWidth', 2, 'Color', colors{i});
    end
    xlabel('Curvature');
    ylabel('Moment');
    title('Moment vs Curvature');
    legend(zone_labels, 'Location', 'best');
    grid on;
    hold off;
end