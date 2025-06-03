function plot_deflection(kappa, epsilon_cr, chi_su, omega, lambda_cu, rho_t, rho_c, beta_all, delta_total, load_4pb, alpha,Envelope, hasFlex,exp_deflection_data, exp_load_data)

    %%%%%%%%%%%%%%%%%%%%%%%%%% Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%
    k_final = Envelope(:,1);
    est_y = kappa * epsilon_cr;
    est_ult = chi_su * epsilon_cr;
    esc_y = omega * epsilon_cr;
    esc_ult = lambda_cu * epsilon_cr;

    % Initialize crossing rows
    crossing_rows = {nan, nan, nan, nan, nan, nan};

    %%%%%%% Steel bottom and steel top yield points %%%%%%%
    est_bot = calculate_est(rho_t, alpha, k_final, beta_all, epsilon_cr, true);
    est_top = calculate_est(rho_c, alpha, k_final, beta_all, epsilon_cr, false);

    %%%%%%% Calculate and Display Yield and Ultimate Points %%%%%%%
    if ~isnan(est_bot)
        [crossing_rows{1}, crossing_rows{2}] = find_crossing_rows(est_bot, est_y, est_ult, beta_all, delta_total, 'Bottom steel');
    end
    if ~isnan(est_top)
        [crossing_rows{3}, crossing_rows{4}] = find_crossing_rows(est_top, est_y, est_ult, beta_all, delta_total, 'Top steel');
    end

    %%%%%%% Concrete yield and ultimate points %%%%%%%
    ectop = k_final .* beta_all * epsilon_cr ./ (1 - k_final);
    [crossing_rows{5}, crossing_rows{6}] = find_crossing_rows(ectop, esc_y, esc_ult, beta_all, delta_total, 'Concrete in compression');

    %%%%%%% Plot the deflection %%%%%%%
    figure;
    if hasFlex == 1
    plot(exp_deflection_data, exp_load_data, '-x', 'DisplayName', 'Experimental Data', 'LineWidth', 1, 'Color', 'r')
    hold on;
    end
    hold on;
    plot(delta_total, load_4pb, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Simulation'); hold on;

    markers = {'ro', 'bo', 'go', 'mo', 'ko', 'co'};
    labels = {'Bottom Steel Yield', 'Bottom Steel Ultimate', 'Top Steel Yield', 'Top Steel Ultimate', 'Concrete Yield', 'Concrete Ultimate'};

    for i = 1:length(crossing_rows)
        if ~isempty(crossing_rows{i})
            % Remove NaN values from rows{i}
            validIndices = crossing_rows{i}(~isnan(crossing_rows{i}));
            if ~isempty(validIndices)
                plot(delta_total(validIndices), load_4pb(validIndices), markers{i}, 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', labels{i});
            else
                % Plot a placeholder point outside the usual data range to ensure it's in the legend
                plot(NaN, NaN, markers{i}, 'MarkerSize', 8, 'LineWidth', 0.2, 'DisplayName', [labels{i}, ' (not found)']);
            end
        else
            % Plot a placeholder point outside the usual data range to ensure it's in the legend
            plot(NaN, NaN, markers{i}, 'MarkerSize', 8, 'LineWidth', 0.2, 'DisplayName', [labels{i}, ' (not found)']);
        end
    end

    legend('Location', 'best');
    xlabel('Deflection');
    ylabel('Load');
    title('Load-Deflection Plot');
    grid on;
    hold on;

end

function est = calculate_est(rho, alpha, k_final, beta_all, epsilon_cr, is_bottom)
    if rho > 0
        if is_bottom
            est = (-alpha + k_final) .* beta_all .* epsilon_cr ./ (k_final - 1);
        else
            est = (k_final - 1 + alpha) .* beta_all * epsilon_cr ./ (k_final - 1);
        end
    else
        est = NaN;
    end
end

function [crossing_row1, crossing_row2] = find_crossing_rows(est, est_y, est_ult, beta_all, delta_total, label)
    crossing_row1 = find(est >= est_y, 1, 'first');
    crossing_row2 = find(est >= est_ult, 1, 'first');

    if ~isempty(crossing_row1)
        beta_yield = beta_all(crossing_row1, 1);
        defl_yield = delta_total(crossing_row1, 1);
        disp([label ' yields at beta: ', num2str(beta_yield)]);
    else
        disp(['No ' label ' yielding point found']);
    end
    if ~isempty(crossing_row2)
        beta_ult = beta_all(crossing_row2, 1);
        defl_ult = delta_total(crossing_row2, 1);
        disp([label ' reaches ultimate strain at beta: ', num2str(beta_ult)]);
    else
        disp(['No ' label ' ultimate point found']);
    end
end