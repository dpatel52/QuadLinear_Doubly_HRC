function [force_Tens,force_Comp, tot_force_rebar_bot, tot_force_rebar_top] = plot_stress_distribution(kd_g, bt_g, d,b,s_area_bot,s_area_top, beta, beta_tu, epsilon_cr, stress, strain, alpha,rho, zeta, strain_st, stress_st)

    n = 5000; %points across the depth
    windowSize = 10;  % Define the window size for the moving average


% Find the index for beta in bt_g, assuming we need to match the closest value
    [~, idx] = min(abs(bt_g - beta));  % Finds the index of the closest match

    % Calculate neutral axis depth using the corresponding kd_g value from the top of the beam
    neutral_axis = kd_g(idx) * d;

    % Strain gradient from neutral axis to the top and bottom of the beam
    % Strain is zero at neutral axis and varies linearly to the top and bottom
    strain_top = -(beta * epsilon_cr) * ( kd_g(idx) / (1- kd_g(idx)));  % Strain at the top of the beam
    strain_bottom = (beta * epsilon_cr);  % Strain at the bottom of the beam


    %make stress strain in steel to for compression as well
    strain_st = [-1*strain_st(3),-1*strain_st(2),strain_st];
    stress_st = [-1*stress_st(3),-1*stress_st(2),stress_st];



    % Creating an array for the depth of the section from top to bottom
    y = linspace(0, d, n);  % n points across the depth

    % Calculate strain at each depth y
    strains = linspace(strain_top, strain_bottom, n);

    % Interpolate stress values for given strains
    stress_interp = griddedInterpolant(strain, stress, 'linear', 'none');
    stress_at_y = stress_interp(strains);
    stress_at_y(isnan(stress_at_y)) = 0;  % Handle NaN values

        % Smooth the stress data to reduce noise
    stress_at_y = smoothdata(stress_at_y, 'movmean', windowSize);

    if rho > 0 && zeta < 0.00001 %condition if it has only bottom rebar
                % Location of rebar (given as alpha relative to depth d)
                rebar_depth = alpha * d;
                rebar_strain = interp1([0, d], [(strain_top), (strain_bottom)], rebar_depth);

                % Interpolating stress for rebar using rebar strain
                rebar_stress_interp = griddedInterpolant(strain_st, stress_st, 'linear', 'none');
                rebar_stress = rebar_stress_interp(rebar_strain);


                % Plotting
                % Fill the regions(only for plotting)
                % Identify positive stress values and corresponding y values
                positive_indices = find(stress_at_y > 0);
                positive_stresses = [0,stress_at_y(positive_indices),0,0];
                y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                % Identify negative stress values and corresponding y values
                negative_indices = find(stress_at_y < 0);
                negative_stresses = [0,0,stress_at_y(negative_indices),0];
                y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];
                figure;
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot(stress_at_y, y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill(positive_stresses, y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill(negative_stresses, y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Stress (psi or MPa)');
                ylabel('Depth from top (inch or mm)');
                title('Stress Distribution');
                grid on;
                hold on;
                % Plot rebar stress
                plot([0,rebar_stress], [rebar_depth,rebar_depth],'Color', 'k', 'LineWidth', 3);

                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');

                % Add neutral axis line at Y=neutral_axis
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Add zero stress line (vertical line)
                %line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

                % Legend and additional plot settings
                %legend({'Stress Profile', 'Tension Zone', 'Compression Zone'});

%plot Force Distribution

                force_at_y = stress_at_y.*b;
                positive_forces = positive_stresses.*b;
                negative_forces = negative_stresses.*b;
                rebar_force_bot =  rebar_stress.*s_area_bot; 

                figure (26) %plot force distribution
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot((force_at_y), y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill((positive_forces), y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill((negative_forces), y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Force (lbf or N)');
                ylabel('Depth from top (inch or mm)');
                title('Force Distribution');
                grid on;
            
                % Plot rebar force
                plot([0, rebar_force_bot], [rebar_depth,rebar_depth],'Color', 'k', 'LineWidth', 3);
                hold on;
                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');
                 hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);


                 %FOR GETTING FORCE DISTRIBUTION
                 for i = 1:length(bt_g)
                        beta = bt_g(i,1)-0.1;

                        [~, idx] = min(abs(bt_g - beta));  % Finds the index of the closest match

                        neutral_axis = kd_g(idx) * d;

                        strain_top = -(beta * epsilon_cr) * ( kd_g(idx) / (1- kd_g(idx)));  
                        strain_bottom = (beta * epsilon_cr); 
                         % Location of rebar (given as alpha relative to depth d)
                         rebar_depth = alpha * d;
                          rebar_strain = interp1([0, d], [(strain_top), (strain_bottom)], rebar_depth);

                       % Interpolating stress for rebar using rebar strain
                      rebar_stress_interp = griddedInterpolant(strain_st, stress_st, 'linear', 'none');
                      rebar_stress = rebar_stress_interp(rebar_strain);

                        %n = 5000; %points across the depth
                        % Creating an array for the depth of the section from top to bottom
                        y = linspace(0, d, n);  % n points across the depth

                        % Calculate strain at each depth y
                        strains = linspace(strain_top, strain_bottom, n);

                        % Interpolate stress values for given strains
                        stress_interp = griddedInterpolant(strain, stress, 'linear', 'none');
                        stress_at_y = stress_interp(strains);
                        stress_at_y(isnan(stress_at_y)) = 0;  % Handle NaN values

        % Smooth stress data
        stress_at_y = smoothdata(stress_at_y, 'movmean', windowSize);

                            positive_indices = find(stress_at_y > 0);
                            positive_stresses = [0,stress_at_y(positive_indices),0,0];
                            y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                            % Identify negative stress values and corresponding y values
                            negative_indices = find(stress_at_y < 0);
                            negative_stresses = [0,0,stress_at_y(negative_indices),0];
                            y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];

                            positive_forces = positive_stresses.*b;
                            negative_forces = negative_stresses.*b;
                            rebar_force_bot =  rebar_stress.*s_area_bot; 
                           % rebar_force_top =  rebar_stress2.*s_area_top;

                           force_Tens(i,1) = trapz(-1*positive_forces,y_positive);
                           force_Comp(i,1) = trapz(-1*negative_forces,y_negative);
                           tot_force_rebar_bot(i,1) = rebar_force_bot;
                           tot_force_rebar_top(i,1) = 0;
                 end


    elseif rho > 0 && zeta > 0.00001 %condition if it has top rebar as well

                % Location of bottom rebar (given as alpha relative to depth d)
                rebar_depth = alpha * d;
                rebar_strain = interp1([0, d], [strain_top, strain_bottom], rebar_depth);

                % Interpolating stress for rebar using rebar strain
                rebar_stress_interp = griddedInterpolant(strain_st, stress_st, 'linear', 'none');
                rebar_stress = rebar_stress_interp(rebar_strain);

                % Location of top rebar (given as alpha relative to depth d)
                rebar_depth2 = (1-alpha) * d;
                rebar_strain2 = interp1([0, d], [strain_top, strain_bottom], rebar_depth2);

                % Interpolating stress for rebar using rebar strain
                rebar_stress_interp2 = griddedInterpolant((strain_st), (stress_st), 'linear', 'none');
                rebar_stress2 = rebar_stress_interp2(rebar_strain2);

                % Plotting
                % Fill the regions(only for plotting)
                % Identify positive stress values and corresponding y values
                positive_indices = find(stress_at_y > 0);
                positive_stresses = [0,stress_at_y(positive_indices),0,0];
                y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                % Identify negative stress values and corresponding y values
                negative_indices = find(stress_at_y < 0);
                negative_stresses = [0,0,stress_at_y(negative_indices),0];
                y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];
                figure;
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot(stress_at_y, y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill(positive_stresses, y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill(negative_stresses, y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Stress (psi or MPa)');
                ylabel('Depth from top (inch or mm)');
                title('Stress Distribution');
                grid on;
                hold on;
                %line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Plot rebar stress
                plot([0,rebar_stress], [rebar_depth,rebar_depth],'Color', 'k', 'LineWidth', 3);
                hold on;
                % Plot rebar stress
                plot([0,rebar_stress2], [rebar_depth2,rebar_depth2],'Color', 'k', 'LineWidth', 3);

                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');
                % Add neutral axis line at Y=neutral_axis
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Add zero stress line (vertical line)
                %line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

                % Legend and additional plot settings
                 %legend({'Stress Profile', 'Tension Zone', 'Compression Zone'});

%plot Force Distribution

                force_at_y = stress_at_y.*b;
                positive_forces = positive_stresses.*b;
                negative_forces = negative_stresses.*b;
                rebar_force_bot =  rebar_stress.*s_area_bot; 
                rebar_force_top =  rebar_stress2.*s_area_top; 

                figure (26) %plot force distribution
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot((force_at_y), y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill((positive_forces), y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill((negative_forces), y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Force (lbf or N)');
                ylabel('Depth from top (inch or mm)');
                title('Force Distribution');
                grid on;
               
                %line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Plot rebar force
                plot([0, rebar_force_bot], [rebar_depth,rebar_depth],'Color', 'k', 'LineWidth', 3);
                hold on;
                % Plot rebar force
                plot([0,rebar_force_top], [rebar_depth2,rebar_depth2],'Color', 'k', 'LineWidth', 3);
                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);


                 %FOR GETTING FORCE DISTRIBUTION
                 for i = 1:length(bt_g)
                        beta = bt_g(i,1)-0.1;

                        [~, idx] = min(abs(bt_g - beta));  % Finds the index of the closest match

                        neutral_axis = kd_g(idx) * d;

                        strain_top = -(beta * epsilon_cr) * ( kd_g(idx) / (1- kd_g(idx)));  
                        strain_bottom = (beta * epsilon_cr); 
                       % Location of rebar (given as alpha relative to depth d)
                         rebar_depth = alpha * d;
                          rebar_strain = interp1([0, d], [(strain_top), (strain_bottom)], rebar_depth);

                       % Interpolating stress for rebar using rebar strain
                      rebar_stress_interp = griddedInterpolant(strain_st, stress_st, 'linear', 'none');
                      rebar_stress = rebar_stress_interp(rebar_strain);
                                      % Location of top rebar (given as alpha relative to depth d)
                rebar_depth2 = (1-alpha) * d;
                rebar_strain2 = interp1([0, d], [strain_top, strain_bottom], rebar_depth2);

                % Interpolating stress for rebar using rebar strain
                rebar_stress_interp2 = griddedInterpolant((strain_st), (stress_st), 'linear', 'none');
                rebar_stress2 = rebar_stress_interp2(rebar_strain2);


                        %n = 5000; %points across the depth
                        % Creating an array for the depth of the section from top to bottom
                        y = linspace(0, d, n);  % n points across the depth

                        % Calculate strain at each depth y
                        strains = linspace(strain_top, strain_bottom, n);

                        % Interpolate stress values for given strains
                        stress_interp = griddedInterpolant(strain, stress, 'linear', 'none');
                        stress_at_y = stress_interp(strains);
                        stress_at_y(isnan(stress_at_y)) = 0;  % Handle NaN values
                                % Smooth stress data
        stress_at_y = smoothdata(stress_at_y, 'movmean', windowSize);


                            positive_indices = find(stress_at_y > 0);
                            positive_stresses = [0,stress_at_y(positive_indices),0,0];
                            y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                            % Identify negative stress values and corresponding y values
                            negative_indices = find(stress_at_y < 0);
                            negative_stresses = [0,0,stress_at_y(negative_indices),0];
                            y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];

                            positive_forces = positive_stresses.*b;
                            negative_forces = negative_stresses.*b;
                            rebar_force_bot =  rebar_stress.*s_area_bot; 
                            rebar_force_top =  rebar_stress2.*s_area_top;

                           force_Tens(i,1) = trapz(-1*positive_forces,y_positive);
                           force_Comp(i,1) = trapz(-1*negative_forces,y_negative);
                           tot_force_rebar_bot(i,1) = rebar_force_bot;
                           tot_force_rebar_top(i,1) =  rebar_force_top;
                 end


    else
%plot Stress Distribtuion
                % Plotting
                % Fill the regions (only for plotting)
                % Identify positive stress values and corresponding y values
                positive_indices = find(stress_at_y > 0);
                positive_stresses = [0,stress_at_y(positive_indices),0,0];
                y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                % Identify negative stress values and corresponding y values
                negative_indices = find(stress_at_y < 0);
                negative_stresses = [0,0,stress_at_y(negative_indices),0];
                y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];
                
                figure;
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot(stress_at_y, y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill(positive_stresses, y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill(negative_stresses, y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Stress (psi or MPa)');
                ylabel('Depth from top (inch or mm)');
                title('Stress Distribution');
                grid on;

                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');

                % Add neutral axis line at Y=neutral_axis
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Add zero stress line (vertical line)
                %line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

                % Legend and additional plot settings
                 %legend({'Stress Profile', 'Tension Zone', 'Compression Zone'});


%plot Force Distribution

                force_at_y = stress_at_y.*b;
                positive_forces = positive_stresses.*b;
                negative_forces = negative_stresses.*b;

                figure (26) %plot force distribution
                subplot(1, 2, 1); % Left subplot for strain distribution
                plot(strains, y, 'Color', 'k', 'LineWidth', 1.5); % Plot strain distribution
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                xlabel('Strain');
                ylabel('Depth from top (inch or mm)');
                title('Strain Distribution');
                line([0 0], ylim, 'Color', 'b', 'LineWidth', 0.5);
                grid on;
                hold on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                set(gca, 'YDir', 'reverse'); % Invert Y-axis so 0 is at the top

                subplot(1, 2, 2); % Right subplot for stress distribution
                plot((force_at_y), y, 'LineWidth', 1); % Thicker line for better visibility
                set(gca, 'YLim', [0 d]); % Set Y-axis limits from 0 to d
                hold on;
                fill((positive_forces), y_positive, 'r', 'FaceAlpha', 0.2); % Filled area with transparency
                hold on;
                fill((negative_forces), y_negative, 'b', 'FaceAlpha', 0.2); % Filled area with transparency
                xlabel('Force (lbf or N)');
                ylabel('Depth from top (inch or mm)');
                title('Force Distribution');
                grid on;
                line(xlim, [neutral_axis neutral_axis], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
                % Invert the Y-axis so 0 is at the top
                set(gca, 'YDir', 'reverse');


                 %FOR GETTING FORCE DISTRIBUTION
                 for i = 1:length(bt_g)
                        beta = bt_g(i,1)-0.1;

                        [~, idx] = min(abs(bt_g - beta));  % Finds the index of the closest match

                        neutral_axis = kd_g(idx) * d;

                        strain_top = -(beta * epsilon_cr) * ( kd_g(idx) / (1- kd_g(idx)));  
                        strain_bottom = (beta * epsilon_cr); 


                        %n = 5000; %points across the depth
                        % Creating an array for the depth of the section from top to bottom
                        y = linspace(0, d, n);  % n points across the depth

                        % Calculate strain at each depth y
                        strains = linspace(strain_top, strain_bottom, n);

                        % Interpolate stress values for given strains
                        stress_interp = griddedInterpolant(strain, stress, 'linear', 'none');
                        stress_at_y = stress_interp(strains);
                        stress_at_y(isnan(stress_at_y)) = 0;  % Handle NaN values
        % Smooth stress data
        stress_at_y = smoothdata(stress_at_y, 'movmean', windowSize);

                            positive_indices = find(stress_at_y > 0);
                            positive_stresses = [0,stress_at_y(positive_indices),0,0];
                            y_positive = [neutral_axis,y(positive_indices),d,neutral_axis];
                            % Identify negative stress values and corresponding y values
                            negative_indices = find(stress_at_y < 0);
                            negative_stresses = [0,0,stress_at_y(negative_indices),0];
                            y_negative = [neutral_axis,0,y(negative_indices),neutral_axis];

                            positive_forces = positive_stresses.*b;
                            negative_forces = negative_stresses.*b;

                           force_Tens(i,1) = trapz(-1*positive_forces,y_positive);
                           force_Comp(i,1) = trapz(-1*negative_forces,y_negative);
                           tot_force_rebar_bot(i,1) = 0;
                           tot_force_rebar_top(i,1) = 0;
                 end


    end
end
