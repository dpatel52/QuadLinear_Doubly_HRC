function draw3DBeam(height, width, length, cover, topDiameter, topCount, botDiameter, botCount, mu, loadingType, loadSpacing)
%%% EXAMPLE (paste this in your program)
% draw3DBeam(h, b, L, cover, topDiameter, topCount, botDiameter, botCount, mu, pointBend, S2); %using US Custom units in recommended

    % Create figure
    figure;
    hold on;
    axis equal;
    view(3); % Set to 3D view

    % Set colors for materials
    concreteColor = [0.7 0.7 0.7]; % Gray
    rebarColor = [0.8 0.2 0.2];    % Reddish
    fiberColor = [0 0 0];          % Black for fibers

    % Draw beam as a transparent block
    [X, Z, Y] = meshgrid([0 length], [0 height], [0 width]);
    V = permute(cat(4, X, Y, Z), [2 1 3 4]);
    patch('Vertices', reshape(V, [], 3), 'Faces', [1 2 6 5; 2 4 8 6; 4 3 7 8; 3 1 5 7; 5 6 8 7; 1 3 4 2], ...
          'FaceColor', concreteColor, 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    % Calculate rebar positions along the width (Y-axis)
      if topCount == 1
        topSpacing = width / 2;
    else
    topSpacing = linspace(cover + topDiameter/2, width - cover - topDiameter/2, topCount);
      end
       if botCount == 1
        botSpacing = width / 2;
    else
    botSpacing = linspace(cover + botDiameter/2, width - cover - botDiameter/2, botCount);
       end
    % Function to draw cylindrical rebars
    function drawRebar(diameter, countSpacing, zPosition)
        for y = countSpacing
            [Xcyl, Ycyl, Zcyl] = cylinder(diameter / 2, 20);
            Zcyl = Zcyl * length;
            Xcyl = Xcyl + zPosition;
            Ycyl = Ycyl + y;
            surf(Zcyl, Ycyl, Xcyl, 'FaceColor', rebarColor, 'EdgeColor', 'none');
            hold on;
            fill3(Zcyl(1, :), Ycyl(1, :), Xcyl(1, :), rebarColor);
            fill3(Zcyl(2, :), Ycyl(2, :), Xcyl(2, :), rebarColor);
        end
    end

    % Draw top and bottom rebars
    drawRebar(topDiameter, topSpacing, height - cover - topDiameter/2);
    drawRebar(botDiameter, botSpacing, cover + botDiameter/2);

    % Add fibers if FRC is enabled

        numFibers = mu*(width*10); % Number of fibers
        fiberLength = 15; % Length of each fiber
        fiberWidth = 0.01; % Thickness of fibers

        for j = 1:numFibers
            valid = false;
            while ~valid
                xStart = rand() * length;
                yStart = rand() * width;
                zStart = rand() * height;
                direction = randn(1,3);
                direction = direction / norm(direction); 
                xEnd = xStart + direction(1) * fiberLength;
                yEnd = yStart + direction(2) * fiberLength;
                zEnd = zStart + direction(3) * fiberLength;
                valid = xEnd >= 0 && xEnd <= length && yEnd >= 0 && yEnd <= width && zEnd >= 0 && zEnd <= height;
            end
            plot3([xStart, xEnd], [yStart, yEnd], [zStart, zEnd], 'Color', fiberColor, 'LineWidth', fiberWidth);
        end


    % Point loads
    floatHeight = height / 2;
    loadArrowLength = height / 2;

    if loadingType == 3
        % Single point load at the center
        quiver3(length/2, width/2, height + floatHeight, 0, 0, -loadArrowLength, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 3, 'AutoScale', 'off');
    else
        % Two point loads with specified spacing
        firstLoadPosition = (length - loadSpacing) / 2;
        secondLoadPosition = firstLoadPosition + loadSpacing;
        quiver3(firstLoadPosition, width/2, height + floatHeight, 0, 0, -loadArrowLength, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 3, 'AutoScale', 'off');
        quiver3(secondLoadPosition, width/2, height + floatHeight, 0, 0, -loadArrowLength, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 3, 'AutoScale', 'off');
    end

    % Additional arrows at the ends and bottom pointing upwards
    arrowHeight = height / 2;
    quiver3(0, width/2, -arrowHeight, 0, 0, arrowHeight, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 3, 'AutoScale', 'off');
    quiver3(length, width/2, -arrowHeight, 0, 0, arrowHeight, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 3, 'AutoScale', 'off');

    % Set plot limits and labels
    xlim([-0.1*length, 1.1*length]);
    ylim([-0.5*width, 1.5*width]);
    zlim([-0.5*height, 1.5*height]);
    xlabel('Length (in)');
    ylabel('Width (in)');
    zlabel('Height (in)');
    title('3D Beam Visualization');
    %grid on;
    hold off;
end
