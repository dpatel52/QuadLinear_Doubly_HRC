function plotBetaVsMandK(beta_all, Envelope, beta_z1, beta_z2, beta_z3, beta_z4, ...
    beta111, M111, k111, beta211, M211, k211, beta212, M212, k212, ...
    beta221, M221, k221, beta222, M222, k222, ...
    beta311, M311, k311, beta312, M312, k312, ...
    beta321, M321, k321, beta322, M322, k322, ...
    beta411, M411, k411, beta412, M412, k412, ...
    beta421, M421, k421, beta422, M422, k422, ...
    beta4222, M4222, k4222)

    % Determine marker spacing
    markerSpacing = size(beta_all, 1) * 0.1; % Change this value to adjust spacing

    % Plot for Beta vs M
    figure;
    hold on;

    % Check and plot each M variable
    if exist('M111', 'var')
        plot(beta111, M111, '-o', 'DisplayName', 'zone111', 'LineWidth', 1, 'Color', 'k', 'MarkerIndices', 1:markerSpacing:length(M111));
    end

    if exist('M211', 'var')
        plot(beta211, M211, '-o', 'DisplayName', 'zone211', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(M211));
    end

    if exist('M212', 'var')
        plot(beta212, M212, '-^', 'DisplayName', 'zone212', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(M212));
    end

    if exist('M221', 'var')
        plot(beta221, M221, '-s', 'DisplayName', 'zone221', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(M221));
    end

    if exist('M222', 'var')
        plot(beta222, M222, '-d', 'DisplayName', 'zone222', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(M222));
    end

    if exist('M311', 'var')
        plot(beta311, M311, '-o', 'DisplayName', 'zone311', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(M311));
    end

    if exist('M312', 'var')
        plot(beta312, M312, '-^', 'DisplayName', 'zone312', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(M312));
    end

    if exist('M321', 'var')
        plot(beta321, M321, '-s', 'DisplayName', 'zone321', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(M321));
    end

    if exist('M322', 'var')
        plot(beta322, M322, '-d', 'DisplayName', 'zone322', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(M322));
    end

    if exist('M411', 'var')
        plot(beta411, M411, '-o', 'DisplayName', 'zone411', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(M411));
    end

    if exist('M412', 'var')
        plot(beta412, M412, '-^', 'DisplayName', 'zone412', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(M412));
    end

    if exist('M421', 'var')
        plot(beta421, M421, '-s', 'DisplayName', 'zone421', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(M421));
    end

    if exist('M422', 'var')
        plot(beta422, M422, '-d', 'DisplayName', 'zone422', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(M422));
    end

    if exist('M4222', 'var')
        plot(beta4222, M4222, '-*', 'DisplayName', 'zone4222', 'LineWidth', 1, 'Color', 'g', 'MarkerIndices', 1:markerSpacing:length(M4222));
    end

    % Add the envelope plot
    beta_plot = [beta_z1; beta_z2; beta_z3; beta_z4];
    plot(beta_plot, Envelope(1:length(beta_plot), 2), 'DisplayName', 'Envelope', 'LineWidth', 5, 'Color', [0 0 0 0.5]);
    ylimit=max(Envelope(1:length(beta_plot), 2));
    hold off;
    legend show;
    xlabel('Beta');
    ylabel('M');
    ylim([0 (ylimit*1.2)]);
    title('Beta vs M for Different Zones');

    % Plot for Beta vs k
    figure;
    hold on;

    % Check and plot each k variable
    if exist('k111', 'var')
        plot(beta111, k111, '-o', 'DisplayName', 'zone111', 'LineWidth', 1, 'Color', 'k', 'MarkerIndices', 1:markerSpacing:length(k111));
    end

    if exist('k211', 'var')
        plot(beta211, k211, '-o', 'DisplayName', 'zone211', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(k211));
    end

    if exist('k212', 'var')
        plot(beta212, k212, '-^', 'DisplayName', 'zone212', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(k212));
    end

    if exist('k221', 'var')
        plot(beta221, k221, '-s', 'DisplayName', 'zone221', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(k221));
    end

    if exist('k222', 'var')
        plot(beta222, k222, '-d', 'DisplayName', 'zone222', 'LineWidth', 1, 'Color', 'm', 'MarkerIndices', 1:markerSpacing:length(k222));
    end

    if exist('k311', 'var')
        plot(beta311, k311, '-o', 'DisplayName', 'zone311', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(k311));
    end

    if exist('k312', 'var')
        plot(beta312, k312, '-^', 'DisplayName', 'zone312', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(k312));
    end

    if exist('k321', 'var')
        plot(beta321, k321, '-s', 'DisplayName', 'zone321', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(k321));
    end

    if exist('k322', 'var')
        plot(beta322, k322, '-d', 'DisplayName', 'zone322', 'LineWidth', 1, 'Color', 'b', 'MarkerIndices', 1:markerSpacing:length(k322));
    end

    if exist('k411', 'var')
        plot(beta411, k411, '-o', 'DisplayName', 'zone411', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(k411));
    end

    if exist('k412', 'var')
        plot(beta412, k412, '-^', 'DisplayName', 'zone412', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(k412));
    end

    if exist('k421', 'var')
        plot(beta421, k421, '-s', 'DisplayName', 'zone421', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(k421));
    end

    if exist('k422', 'var')
        plot(beta422, k422, '-d', 'DisplayName', 'zone422', 'LineWidth', 1, 'Color', 'r', 'MarkerIndices', 1:markerSpacing:length(k422));
    end

    if exist('k4222', 'var')
        plot(beta4222, k4222, '-*', 'DisplayName', 'zone4222', 'LineWidth', 1, 'Color', 'g', 'MarkerIndices', 1:markerSpacing:length(k4222));
    end

    hold off;
    xlabel('Beta');
    ylabel('k');
    title('Beta vs k for Different Zones');
    legend show;
end
