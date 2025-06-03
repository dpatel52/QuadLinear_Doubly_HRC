function drawDoublyReinforcedBeam(height, width, cover, topDiameter, topCount, botDiameter, botCount,mu)

%%%%% Example (paste this in your code)
%drawDoublyReinforcedBeam(h, b, cover, topDiameter, topCount, botDiameter, botCount,mu); %using US Custom units in recommended

    % Create figure
    figure;
    hold on;
    axis equal;
    %grid on;

    % Set colors for materials
    concreteColor = [0.7 0.7 0.7]; % Gray
    rebarColor = [0.8 0.2 0.2];    % Reddish
    fiberColor = [0 0 0];          % Black for fibers

    % Draw concrete cross section
    rectangle('Position', [0, 0, width, height], 'FaceColor', concreteColor, 'EdgeColor', 'k');

    % Calculate rebar positions
    topRebarY = height - cover;
    botRebarY = cover;

  

    % % if topCount > 1
    % %     topSpacing = (width - 2*cover - topDiameter) / (topCount - 1);
    % % else
    % %     topSpacing = 0;
    % % end
    % % 
    % % if botCount > 1
    % %     botSpacing = (width - 2*cover - botDiameter) / (botCount - 1);
    % % else
    % %     botSpacing = 0;
    % % end

  

    % If FRC is set, draw random fibers
        % Define number of fibers and their length
        if width >74 || height > 74
        numFibers = mu*300; % Adjust number of fibers as needed
        else
        numFibers = mu*(height*width)*1.1; % Adjust number of fibers as needed
        end
        if width >74 || height > 74
        fiberLength = 10; % Fixed length mm ( 1 inch )
        else
        fiberLength = 1; % Fixed length of 1 inch
        end
        fiberWidth = 0.01; % Reducing fiber width for a more realistic appearance

        for j = 1:numFibers
            % Generate a random angle in radians
            angle = rand(1) * 2 * pi; % Fibers can be oriented in any direction in 360 degrees

            % Calculate the fiber endpoints based on angle and length
            deltaX = cos(angle) * fiberLength;
            deltaY = sin(angle) * fiberLength;

            % Define the allowable starting coordinates to ensure the fiber stays within boundaries
            xMin = max(0, -deltaX); % If deltaX is negative, adjust xMin accordingly
            yMin = max(0, -deltaY); % If deltaY is negative, adjust yMin accordingly

            xMax = min(width, width - deltaX); % Adjust xMax based on deltaX to keep within width
            yMax = min(height, height - deltaY); % Adjust yMax based on deltaY to keep within height

            % Random starting position ensuring both ends of the fiber are inside the beam
            xStart = xMin + (xMax - xMin) * rand(1);
            yStart = yMin + (yMax - yMin) * rand(1);

            % Calculate the ending position
            xEnd = xStart + deltaX;
            yEnd = yStart + deltaY;

            % Plot each fiber
            plot([xStart, xEnd], [yStart, yEnd], 'Color', fiberColor, 'LineWidth', fiberWidth);
        end
          % Draw top rebars
    % for i = 0:topCount-1
    %     x = cover + topDiameter/2 + i*topSpacing;
    %     rectangle('Position', [x - topDiameter/2, topRebarY - topDiameter/2, topDiameter, topDiameter], ...
    %               'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
    % end
    % 
    % % Draw bottom rebars
    % for i = 0:botCount-1
    %     x = cover + botDiameter/2 + i*botSpacing;
    %     rectangle('Position', [x - botDiameter/2, botRebarY - botDiameter/2, botDiameter, botDiameter], ...
    %               'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
    % end
       
  % Even spacing for top and bottom rebars
    if topCount > 1
        topSpacing = (width - 2 * (cover + topDiameter / 2)) / (topCount - 1);
        for i = 0:topCount-1
            x = cover + topDiameter / 2 + i * topSpacing;
            rectangle('Position', [x - topDiameter / 2, topRebarY - topDiameter / 2, topDiameter, topDiameter], ...
                      'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
        end
    elseif topCount == 1
        % Center the single top rebar
        x = width / 2;
        rectangle('Position', [x - topDiameter / 2, topRebarY - topDiameter / 2, topDiameter, topDiameter], ...
                  'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
    else
        topSpacing = 0;
         for i = 0:topCount-1
        x = cover + topDiameter/2 + i*topSpacing;
        rectangle('Position', [x - topDiameter/2, topRebarY - topDiameter/2, topDiameter, topDiameter], ...
                  'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
         end
    end

    if botCount > 1
        botSpacing = (width - 2 * (cover + botDiameter / 2)) / (botCount - 1);
        for i = 0:botCount-1
            x = cover + botDiameter / 2 + i * botSpacing;
            rectangle('Position', [x - botDiameter / 2, botRebarY - botDiameter / 2, botDiameter, botDiameter], ...
                      'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
        end
    elseif botCount == 1
          % Center the single bottom rebar
        x = width / 2;
        rectangle('Position', [x - botDiameter / 2, botRebarY - botDiameter / 2, botDiameter, botDiameter], ...
                  'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
    else
        botSpacing = 0;

        for i = 0:botCount-1
            x = cover + botDiameter/2 + i*botSpacing;
            rectangle('Position', [x - botDiameter/2, botRebarY - botDiameter/2, botDiameter, botDiameter], ...
                      'Curvature', [1, 1], 'FaceColor', rebarColor, 'EdgeColor', 'k');
        end
    end


        % Set plot limits and labels
    xlim([-cover, width + cover]);
    ylim([-cover, height + cover]);
    xlabel('Width (in)');
    ylabel('Height (in)');
    title('Beam Cross-Section');
    hold off;
end