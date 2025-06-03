function draw2DBeam(height, width, length, cover, topDiameter, topCount, botDiameter, botCount, mu, loadingType, loadSpacing)
    % Clear current figure
    clf;
    hold on;
    axis equal;

    % Set colors for materials
    concreteColor = [0.7 0.7 0.7]; % Gray
    rebarColor = [0.8 0.2 0.2];    % Reddish
    fiberColor = [0 0 0];          % Black for fibers

    % Draw beam outline
    rectangle('Position', [0 0 length height], 'FaceColor', concreteColor, 'EdgeColor', 'k');

    % Draw top rebars
    if topCount > 0 && topDiameter > 0
        topY = height - cover - topDiameter;
        rectangle('Position', [0, topY, length, topDiameter], ...
                  'FaceColor', rebarColor, 'EdgeColor', 'none');
    end

    % Draw bottom rebars
    if botCount > 0 && botDiameter > 0
        botY = cover;
        rectangle('Position', [0, botY, length, botDiameter], ...
                  'FaceColor', rebarColor, 'EdgeColor', 'none');
    end

    % Add fibers
    numFibers = mu * 200;
    fiberLength = height / 10;

    for j = 1:numFibers
        xStart = rand() * length;
        yStart = rand() * height;
        angle = rand() * 2 * pi;
        xEnd = xStart + cos(angle) * fiberLength;
        yEnd = yStart + sin(angle) * fiberLength;

        if xEnd >= 0 && xEnd <= length && yEnd >= 0 && yEnd <= height
            plot([xStart, xEnd], [yStart, yEnd], 'Color', fiberColor, 'LineWidth', 0.5);
        end
    end

    % Draw loads
    arrowLength = height / 4;

    if loadingType == 3
        % Single point load at center
        drawArrow(length/2, height + arrowLength, 0, -arrowLength);
    else
        % Two point loads
        firstLoadPosition = (length - loadSpacing) / 2;
        secondLoadPosition = firstLoadPosition + loadSpacing;
        drawArrow(firstLoadPosition, height + arrowLength, 0, -arrowLength);
        drawArrow(secondLoadPosition, height + arrowLength, 0, -arrowLength);
    end

    % Draw support arrows
    supportArrowLength = height / 4;
    drawArrow(0, -supportArrowLength, 0, supportArrowLength);
    drawArrow(length, -supportArrowLength, 0, supportArrowLength);
    drawArrow(length + supportArrowLength, 0, -supportArrowLength, 0);

    % Set axis limits and grid
    xlim([-0.1*length, 1.2*length]);
    ylim([-0.4*height, 1.5*height]);
    xlabel('Length')
    ylabel('Height')
    grid on;
    hold off;
end

% Helper function to draw arrows
function drawArrow(x, y, dx, dy)
    quiver(x, y, dx, dy, 0, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 1);
end
