clear;
close all;
clc;

% Start profiling
%profile on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       BACK-CALCULATION with QUAD-LINEAR UHPC MODEL W/ DBL. REINF.                                                                                                                                                                                                                                                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb=matlab.codetools.requiredFilesAndProducts('Main_Quad_Doubly.m')';

%Version 9
% Author: Devansh Patel
% Date: June 2025 (Github Version)

% Notes:
% - Moment needs fixing

Create_Outputfile = 0;
Create_pdf = 0;
plot_stress_dist = 0;
output_filename = '********.xlsx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Experimental Data Extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hasFlex = 1; % Set to 1 if experimental data is available, otherwise set to 0
exp_load_data =[];  exp_deflection_data=[];
if hasFlex == 1
    try
        data = readtable('****.xlsx');
        exp_load_data = data{:, 2};      % Load data in the first column
        exp_deflection_data = data{:, 1}; % Deflection data in the second column
    catch
        disp('Error: Could not read the Excel file. Please check the file and try again.');
    end
else
    disp('No experimental data will be plotted.');
end

hasTens = 1; % Set to 1 if experimental data is available, otherwise set to 0
tens_strs_data =[];  tens_strain_data=[];
if hasFlex == 1
    try
        dataT = readtable('****.xlsx');
        tens_strs_data = dataT{:, 2};      % Strs data in the first column
        tens_strain_data = dataT{:, 1}; % Strn data in the second column
    catch
        disp('Error: Could not read the Excel file for Tension Test. Please check the file and try again.');
    end
else
    disp('No experimental Tension data will be plotted.');
end

hasComp = 1; % Set to 1 if experimental data is available, otherwise set to 0
comp_strs_data =[];  comp_strain_data=[];
if hasComp == 1
    try
        dataC = readtable('*****.xlsx');
        comp_strs_data = dataC{:, 2};      % Strs data in the second column
        comp_strain_data = dataC{:, 1}; % Strn data in the first column
    catch
        disp('Error: Could not read the Excel file for Compression Test. Please check the file and try again.');
    end
else
    disp('No experimental Compression data will be plotted.');
end

hasReinf = 1; % Set to 1 if experimental data is available, otherwise set to 0
reinf_strs_data =[];  reinf_strain_data=[];
if hasReinf == 1
    try
        dataR = readtable('****.xlsx');
        reinf_strs_data = dataR{:, 2};      % Strs data in the second column
        reinf_strain_data = dataR{:, 1}; % Strn data in the first column
    catch
        disp('Error: Could not read the Excel file for Reinforcement Test. Please check the file and try again.');
    end
else
    disp('No experimental Reinforcement data will be plotted.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Input Paremeters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================= Geometry and Loading ===========================%

pointBend = 4;        % Number of point loads applied
L = (1500);               % Length of the beam (in inches)
b = 150;                % Width of the beam (in inches)
h = 500;                % Height of the beam (in inches)
S2 = 125;               % Spacing between Loading Points (4PB) only
Lp = 125;               % Length of the span (in inches)
PostLp = 35;          % Half of the height (in inches)
cover = 25;            % Concrete cover (in inches)
alpha = (h - cover) / h; % Ratio of effective height

%==================== Tension Model Parameters ==========================%

fc = 149;
E = 34000;          % Modulus of elasticity (psi)
epsilon_cr = 0.00012; %0.3*(fck)^(2/3)/E;% Cracking strain
mu_1 = 0.33;         % Multiplier for post crack strength at beta 1
mu_2 = 0.4;         % Multiplier for post crack strength at beta 2
mu_3 = 0.33;       % Multiplier for post crack strength at beta 3
beta_1 = 4;          % strain parameter for zone 1
beta_2 = 60;          % strain parameter for zone 2
beta_3 = 150;         % strain parameter for zone 3
eta_1 = (mu_1 - 1) / (beta_1 - 1); % Slope for zone 1
eta_2 = (mu_2 - mu_1) / (beta_2 - beta_1); % Slope for zone 2
eta_3 =(mu_3 - mu_2) / (beta_3 - beta_2); % Slope for zone 3
%mu_3 = eta_3*(beta_3-beta_2)+mu_2;


%=================== Compression Model Parameters =======================%

%Ec = 8010*(fc^0.36);
xi = 1.01;           % Ex/E ratio (>1.01)
%ecy = 0.003;%fc/(E*xi);
omega = 7.3;         % Transition point for compression model
mu_c = 1;           % Multiplier for compression
ecu = 0.003;
lambda_cu = ecu/epsilon_cr;      % Ultimate strain in compression
eta_c = (mu_c - 1) / (lambda_cu - omega); % Slope for compression model

%======================== Steel Properties ==============================%%

Es = 200000;
n = Es/E;                % Es/E ratio
kappa = 0.002/epsilon_cr;           % Nom. Yield strain esy/ecr
mu_s = 1;           % Multiplier for steel
chi_su = 6 * kappa;   % Ultimate strain in steel
eta_s = (mu_s - 1) / (chi_su - kappa); % Slope for steel model

%======================= Reinforcement Details ==========================%

topDiameter = 12;  % Diameter of top bars (inches) - #3 bar
topCount = 0;         % Number of top bars
botDiameter = 18;  % Diameter of bottom bars (inches) - #4 bar
botCount = 0;         % Number of bottom bars

%====================== Reinforcement Areas =============================%

s_area_bot = botCount * (botDiameter^2 * pi / 4); % Bottom steel area
rho_t = s_area_bot / (b * h);  % Reinforcement ratio in tension
s_area_top = topCount * (topDiameter^2 * pi / 4); % Top steel area
rho_c = s_area_top / (b * h);  % Reinforcement ratio in compression

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Draw Cross-Sections                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make 3D plot 
mu_avg = (mu_1 + mu_2 + mu_3) / 3;
%draw3DBeam(h, b, L, cover, topDiameter, topCount, botDiameter, botCount, mu_avg, pointBend, Lp); % using US Custom units as recommended
draw2DBeam(h, b, L, cover, topDiameter, topCount, botDiameter, botCount, mu_avg, pointBend, S2);
drawDoublyReinforcedBeam(h, b, cover, topDiameter, topCount, botDiameter, botCount, mu_avg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Plot Material Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tension Model (Matrix)
strainT = [0, epsilon_cr, epsilon_cr * beta_1, epsilon_cr * beta_2, epsilon_cr * beta_3];
stressT = [0, epsilon_cr * E, mu_1 * epsilon_cr * E, mu_2 * epsilon_cr * E, mu_3 * epsilon_cr * E];

% Compression Model (Matrix)
strainC = [0, omega * epsilon_cr, epsilon_cr * lambda_cu];
stressC = [0, epsilon_cr * E * xi * omega, epsilon_cr * E * xi * omega * mu_c];

% Reinforcement Model
strainR = [0, kappa * epsilon_cr, epsilon_cr * chi_su];
stressR = [0, epsilon_cr * E * n * kappa, epsilon_cr * E * n * kappa * mu_s];

% Plot Tension Model
figure;
plot(strainT, stressT, '-o', 'LineWidth', 2, 'Color', 'r'); % Red for Tension Model
title('Tension Model (Matrix)','FontSize', 14, 'FontWeight', 'bold');
xlabel('Strain (mm/mm)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Stress (MPa)','FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14);
grid on;

% Plot Compression Model
figure;
plot(strainC, stressC, '-o', 'LineWidth', 2, 'Color', 'b'); % Blue for Compression Model
title('Compression Model (Matrix)','FontSize', 14, 'FontWeight', 'bold');
xlabel('Strain (mm/mm)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Stress (MPa)','FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14);
grid on;

% Plot Reinforcement Model
figure;
plot(strainR, stressR, '-o', 'LineWidth', 2, 'Color', [0 0.5 0]); % Dark Green for Reinforcement Model
title('Reinforcement Model','FontSize', 14, 'FontWeight', 'bold');
xlabel('Strain (mm/mm)','FontSize', 14, 'FontWeight', 'bold');
ylabel('Stress (MPa)','FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14); % Increase font size for numbers on the x and y axis
grid on;
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Cracking Properties                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain k from zone111 function
[~, kcr, ~, ~, ~] = zone111(1, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t); 

% Calculate cracking moment, curvature, and flexural rigidity
M_cr = (epsilon_cr * E * b * h^2) / (12 * (1 - kcr)); 
Phi_cr = epsilon_cr / ((1 - kcr) * h);
EI_cr = M_cr / Phi_cr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Beta Arrays Initialization                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create segments with different point densities
%segment1 = linspace(0, beta_3, 500)';                  % 1000 points from 0 to 1

% Combine all segments
beta_z1 = linspace(0, 1, 500)';
beta_z2 = linspace(1, beta_1, 1000)';
beta_z3 = linspace(beta_1, beta_2, 1000)';
beta_z4 = linspace(beta_2, beta_3, 500)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Zone Analysis and Calculations                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis for different zones
[beta111,k111,M111,phi111,defl111,defl_cant111] = zone111(beta_z1,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta211,k211,M211,phi211,defl211,defl_cant211,ectop211,es_T211,yield_type211] = zone211(beta_z2,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta212,k212,M212,phi212,defl212,defl_cant212,ectop212,es_T212,yield_type212] = zone212(beta_z2,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta221,k221,M221,phi221,defl221,defl_cant221,ectop221,es_T221,yield_type221] = zone221(beta_z2,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta222,k222,M222,phi222,defl222,defl_cant222,ectop222,es_T222] = zone222(beta_z2,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta311,k311,M311,phi311,defl311,defl_cant311,ectop311,es_T311,yield_type311] = zone311(beta_z3,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta312,k312,M312,phi312,defl312,defl_cant312,ectop312,es_T312,yield_type312] = zone312(beta_z3,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta321,k321,M321,phi321,defl321,defl_cant321,ectop321,es_T321,yield_type321] = zone321(beta_z3,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta322,k322,M322,phi322,defl322,defl_cant322,ectop322,es_T322] = zone322(beta_z3,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta411,k411,M411,phi411,defl411,defl_cant411,ectop411,es_T411,yield_type411] = zone411(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta412,k412,M412,phi412,defl412,defl_cant412,ectop412,es_T412,yield_type412] = zone412(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta421,k421,M421,phi421,defl421,defl_cant421,ectop421,es_T421,yield_type421] = zone421(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta422,k422,M422,phi422,defl422,defl_cant422,ectop422,es_T422,es_C422,yield_type422] = zone422(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);
[beta4222,k4222,M4222,phi4222,defl4222,defl_cant4222,ectop4222,es_T4222,es_C4222] = zone4222(beta_z4,L,b,h,alpha,E,epsilon_cr,beta_1,beta_2,beta_3,eta_1,eta_2,eta_3,xi,omega,eta_c,n,kappa,eta_s,rho_c,rho_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Yield Point Analysis and Calculations                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_all=[beta_z1;beta_z2;beta_z3;beta_z4];

Envelope = Envelope_Final (kappa, omega, epsilon_cr,beta_z1,beta_z2,beta_z3,beta_z4, k111, M111, ...
                                              k211, M211, k212, M212, k221, M221, k222, M222, ...
                                              k311, M311, k312, M312, k321, M321, k322, M322, ...
                                              k411, M411, k412, M412, k421, M421, k422, M422, ...
                                              k4222,M4222,beta_1,beta_2,beta_3,alpha);
Envelope_beta=[beta_all,Envelope];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Moment & k Plotter at Diffrent Zones                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 plotBetaVsMandK(beta_all, Envelope, beta_z1, beta_z2, beta_z3, beta_z4, ...
    beta111, M111, k111, beta211, M211, k211, beta212, M212, k212, ...
    beta221, M221, k221, beta222, M222, k222, ...
    beta311, M311, k311, beta312, M312, k312, ...
    beta321, M321, k321, beta322, M322, k322, ...
    beta411, M411, k411, beta412, M412, k412, ...
    beta421, M421, k421, beta422, M422, k422, ...
    beta4222, M4222, k4222);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Moment vs Curvature Plotter                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strain_bot = beta_all*epsilon_cr;
NA_from_bot = (1-Envelope(:,1))*h;
Curvature = strain_bot./NA_from_bot;
plot_moment_vs_curvature(Curvature, Envelope, beta_z1, beta_z2, beta_z3, beta_z4)

Moment=Envelope(:,2);
k_final = Envelope(:,1);
est_bot = (-alpha + k_final) .* beta_all .* epsilon_cr ./ (k_final - 1);
% Calculate limit states on the curvature
%service_curv = crackwidth(find(est_bot >= 0.8 * kappa * epsilon_cr, 1), 1);
service_mom = Moment(find(est_bot >= 0.8 * kappa * epsilon_cr, 1), 1);

% Display the values using disp
%disp(['Service Crackwidth: ', num2str(service_curv)]);
disp(['Service Moment: ', num2str(service_mom)]);

k_final = Envelope(:,1);
est_bot = (-alpha + k_final) .* beta_all .* epsilon_cr ./ (k_final - 1);
est_top = abs((k_final - 1 + alpha) .* beta_all * epsilon_cr ./ (k_final - 1));
strain_top = k_final.*beta_all*epsilon_cr./(1 - k_final);

% Define whether the plot is in kip-ft or lb-in
plot_in_kip_ft = 1; % Set to 0 for lb-in, 1 for kip-ft or N-mm to kn-m

% Calculate limit states on the curvature
service_curv = Curvature(find(est_bot >= 0.8 * kappa * epsilon_cr, 1), 1);
bot_st_yld_curv = Curvature(find(est_bot >= kappa * epsilon_cr, 1), 1);
bot_st_ult_curv = Curvature(find(est_bot >= chi_su * epsilon_cr, 1), 1);
top_st_yld_curv = Curvature(find(est_top >= kappa * epsilon_cr, 1), 1);
top_st_ult_curv = Curvature(find(est_top >= chi_su * epsilon_cr, 1), 1);
crack_curv = Phi_cr;
crack_local_curv = Curvature(Envelope(:, 2) == max(Envelope(:, 2)), 1);
comp_yld_curv = Curvature(find(strain_top >= omega * epsilon_cr, 1), 1);
comp_ult_curv = Curvature(find(strain_top >= lambda_cu * epsilon_cr, 1), 1);

% Conversion factor for moment if plotting in kip-ft or kN-m
moment_conversion_factor = 1;
y_label_text = 'Bending Moment (N-mm)';

if plot_in_kip_ft
    moment_conversion_factor = 1 / 12000; % Convert from lb-in to kip-ft
    y_label_text = 'Bending Moment (kip-ft)';
end
hasMomCrv =0;
% Plot Moment vs Curvature
figure;
if hasMomCrv == 1
    plot(curv_data, mom_data * moment_conversion_factor, '-x', 'DisplayName', 'Experimental', 'LineWidth', 1, 'Color', 'r');
end
hold on;
plot(Curvature, (Envelope(:, 2)) * moment_conversion_factor, '-r', 'DisplayName', 'Simulation', 'LineWidth', 2);
hold on;

% Define the markers with better colors for a white background
markers = {'ro', 'bo', 'go', 'mo', 'co', 'ko', 'rs', 'bs', 'ms'};

% Plot markers for each curvature point if it exists, using the specified markers
if ~isempty(crack_curv)
    idx = find(Curvature == crack_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{1}, 'DisplayName', 'First Crack', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(service_curv)
    idx = find(Curvature == service_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{6}, 'DisplayName', 'Service Stress', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(bot_st_yld_curv)
    idx = find(Curvature == bot_st_yld_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{2}, 'DisplayName', 'Steel Yield', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(bot_st_ult_curv)
    idx = find(Curvature == bot_st_ult_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{3}, 'DisplayName', 'Steel Ultimate Strain', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(top_st_yld_curv) && rho_c > 0
    idx = find(Curvature == top_st_yld_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{4}, 'DisplayName', 'Top Steel Yield', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(top_st_ult_curv)
    idx = find(Curvature == top_st_ult_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{5}, 'DisplayName', 'Top Steel Ultimate', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(crack_local_curv)
    idx = find(Curvature == crack_local_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{7}, 'DisplayName', 'Crack Localization', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(comp_yld_curv)
    idx = find(Curvature == comp_yld_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{8}, 'DisplayName', 'Conc. Compression Yield', 'MarkerSize', 12, 'LineWidth', 2);
end
if ~isempty(comp_ult_curv)
    idx = find(Curvature == comp_ult_curv, 1);
    plot(Curvature(idx), Envelope(idx, 2) * moment_conversion_factor, markers{9}, 'DisplayName', 'Conc. Compression Ultimate', 'MarkerSize', 12, 'LineWidth', 2);
end

% Add labels and legend
xlabel('Sectional Curvature (1/in)','FontSize', 14, 'FontWeight', 'bold');
ylabel(y_label_text, 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
%title('Moment vs Curvature');
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Deflection Calculation and Load vs Deflection Plotter              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pointBend == 4
    S1=(L-S2)/2;
load = (Envelope(:,2)*2)/S1;
elseif pointBend == 3
    load = (Envelope(:,2)*4)/L;
else 
    load = (Envelope(:,2)*4)/L;
end
%defl_4pb=-1*(beta_all .* x1 .* epsilon_cr .* (-x1 + L) ./ (2 .* (Envelope(:,1) - 1) .* h));



%Moment-Aread Method for getting deflection
MC_data = [beta_all,Envelope(:,2),Curvature];
Mmax = max(MC_data(:,2)); % Maximum moment at the loading points
idx = find(MC_data(:,2) == Mmax, 1);  % Finds the first occurrence where moment equals Mmax
Cmax = MC_data(idx, 3); %curvature at max moment

mom = MC_data(:,2);
cv=MC_data(:,3);

% Load vs Deflection plot
markerSpacing = size(beta_all, 1) * 0.1;
markerSpacing2 = round(markerSpacing / 5);


%%%%%%%%%%%%%%%% for cantilever

% position_factor = 1;
% bb = position_factor*L;
% delta_total = calculate_deflection_canti(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, PostLp, Mmax,position_factor);
% Load_cant = mom/bb;
% 
% figure;
%     plot(delta_total, Load_cant, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Simulation'); hold on;
%     xlabel('Deflection');
%     ylabel('Load');
%     title('Load-Deflection Plot');
%     grid on;
%     hold on;

%%%%%%%%%%%%%% for 3PB and 4PB
if pointBend == 4
delta_total = calculate_deflection4PB(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, PostLp, Mmax,S2);
elseif pointBend == 3
delta_total = calculate_deflection_3PB(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, PostLp, Mmax,S2);
else 
delta_total = calculate_deflection4PB(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, PostLp, Mmax,S2);
end

%delta_total = calculate_deflection(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, PostLp, Mmax,S2);
plot_deflection(kappa, epsilon_cr, chi_su, omega, lambda_cu, rho_t, rho_c, beta_all, delta_total, load, alpha,Envelope, hasFlex,exp_deflection_data, exp_load_data)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                          Plot Stress Distribution                           
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_stress_dist == 1

stress = [-epsilon_cr*E*xi*omega*mu_c,-epsilon_cr*E*xi*omega,0,epsilon_cr*E,mu_1*epsilon_cr*E, mu_2*epsilon_cr*E,mu_3*epsilon_cr*E];
strain = [-1*(epsilon_cr*lambda_cu),-1*(omega*epsilon_cr), 0, epsilon_cr, epsilon_cr*beta_1, epsilon_cr*beta_2, epsilon_cr*beta_3];

beta_strs = beta_3;
kd_g = Envelope(:,1);
[force_Tens,force_Comp, tot_force_rebar_bot, tot_force_rebar_top] = plot_stress_distribution(kd_g, beta_all, h,b,s_area_bot,s_area_top, beta_strs, beta_3, epsilon_cr, stress, strain, alpha,rho_t, (rho_c/rho_t), strainR, stressR);

beta_final = beta_all;
figure(27); % Open figure window with identifier 26
plot(beta_final, force_Tens, 'r', 'LineWidth', 2); % Plot force_Tens with a blue line
hold on; % Keep the plot for adding more plots or lines
plot(beta_final, force_Comp, 'b', 'LineWidth', 2); % Plot force_Comp with a red line
hold on;
plot(beta_final, tot_force_rebar_bot ,'r', 'LineWidth', 2,'LineStyle', '--'); % Plot force_Comp with a red line
hold on;
plot(beta_final, tot_force_rebar_top, 'b', 'LineWidth', 2,'LineStyle', '--'); % Plot force_Comp with a red line

% % Check if 0 is within the range of beta_final to add a vertical line
% if min(beta_final) <= 0 && max(beta_final) >= 0
%     line([0 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5); % Add vertical line at x = 0
% end

xlabel('Normalized Tensile Strain (beta)'); % Label for the x-axis
ylabel('Force (lbf)'); % Label for the y-axis
%title('Force'); % Title of the plot
legend('Concrete in Tension', 'Concrete in Compression','Rebar in Tension','Rebar in Compression', 'Location', 'best'); % Add a legend
grid on; % Turn on the grid for easier readability
hold on; % Release the figure for other plots


%%%Efficiency

Con_Tens_eff = abs(force_Tens)./(abs(force_Tens)+abs(tot_force_rebar_bot));
Tens_rebar_eff = abs(tot_force_rebar_bot)./(abs(force_Tens)+abs(tot_force_rebar_bot));
Con_Comp_eff = abs(force_Comp)./(abs(force_Comp)+abs(tot_force_rebar_top));
Comp_rebar_eff = abs(tot_force_rebar_top)./(abs(force_Comp)+abs(tot_force_rebar_top));

figure (30);

plot(beta_final, Con_Tens_eff, 'r', 'LineWidth', 2); % Plot force_Tens with a blue line
hold on; % Keep the plot for adding more plots or lines
plot(beta_final, Con_Comp_eff, 'b', 'LineWidth', 2); % Plot force_Comp with a red line
hold on;
plot(beta_final, Tens_rebar_eff ,'r', 'LineWidth', 2,'LineStyle', '--'); % Plot force_Comp with a red line
hold on;
plot(beta_final, Comp_rebar_eff , 'b', 'LineWidth', 2,'LineStyle', '--'); % Plot force_Comp with a red line


xlabel('Normalized Tensile Strain (beta)'); % Label for the x-axis
ylabel('Efficiency factors'); % Label for the y-axis
%title('Force'); % Title of the plot
legend('Concrete in Tension', 'Concrete in Compression','Rebar in Tension','Rebar in Compression', 'Location', 'best'); % Add a legend
grid on; % Turn on the grid for easier readability
hold on; % Release the figure for other plots
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Generate pdf report                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Create_pdf == 1

% Create a new report
import mlreportgen.report.*
import mlreportgen.dom.*

% Prompt the user for a file name
fileName = input('Enter the desired PDF file name (without extension): ', 's');
fileName = strcat(fileName, '.pdf');

% Create the report
rpt = Report(fileName, 'pdf');
open(rpt);

% Add a title page
tp = TitlePage;
tp.Title = 'Back-Calculation with Quad-Linear UHPC Model with Double Reinforcement';
tp.Author = 'Devansh Patel';
tp.Publisher = 'June 2024';
add(rpt, tp);

% Add a single chapter for all content
ch = Chapter('Report');

% Add input parameters section
sec1 = Section('Input Parameters');

% Geometry and Loading table
geoTitle = Paragraph('Geometry and Loading');
geoTitle.Style = {Bold(true)};
add(sec1, geoTitle);

geoTable = {...
    'Number of point loads applied', pointBend; 
    'Length of the beam (in inches)', L; 
    'Width of the beam (in inches)', b;
    'Height of the beam (in inches)', h;
    'Spacing between Loading Points (in inches)', S2;
    'Length of the span (in inches)', Lp;
    'Concrete cover (in inches)', cover;
    'Ratio of effective height', alpha};
t1 = FormalTable(geoTable);
t1.Style = {Border('solid'), RowSep('solid'), ColSep('solid')};
add(sec1, t1);

% Tension Model Parameters table
tensionTitle = Paragraph('Tension Model Parameters');
tensionTitle.Style = {Bold(true)};
add(sec1, tensionTitle);

tensionTable = {...
    'Modulus of elasticity (psi)', E;
    'Cracking strain', epsilon_cr;
    'Multiplier for post crack strength at beta 1', mu_1;
    'Multiplier for post crack strength at beta 2', mu_2;
    'Multiplier for post crack strength at beta 3', mu_3;
    'Strain parameter for zone 1', beta_1;
    'Strain parameter for zone 2', beta_2;
    'Strain parameter for zone 3', beta_3;
    'Slope for zone 1', eta_1;
    'Slope for zone 2', eta_2;
    'Slope for zone 3', eta_3};
t2 = FormalTable(tensionTable);
t2.Style = {Border('solid'), RowSep('solid'), ColSep('solid')};
add(sec1, t2);

% Compression Model Parameters table
compressionTitle = Paragraph('Compression Model Parameters');
compressionTitle.Style = {Bold(true)};
add(sec1, compressionTitle);

compressionTable = {...
    'Compressive strength factor', xi;
    'Transition point for compression model', omega;
    'Multiplier for compression', mu_c;
    'Ultimate strain in compression', lambda_cu;
    'Slope for compression model', eta_c};
t3 = FormalTable(compressionTable);
t3.Style = {Border('solid'), RowSep('solid'), ColSep('solid')};
add(sec1, t3);

% Reinforcement Details table
reinforcementTitle = Paragraph('Reinforcement Details');
reinforcementTitle.Style = {Bold(true)};
add(sec1, reinforcementTitle);

reinforcementTable = {...
    'Number of reinforcing bars', n;
    'Strain hardening factor', kappa;
    'Multiplier for steel', mu_s;
    'Ultimate strain in steel', chi_su;
    'Slope for steel model', eta_s;
    'Diameter of top bars (inches)', topDiameter;
    'Number of top bars', topCount;
    'Diameter of bottom bars (inches)', botDiameter;
    'Number of bottom bars', botCount};
t4 = FormalTable(reinforcementTable);
t4.Style = {Border('solid'), RowSep('solid'), ColSep('solid')};
add(sec1, t4);

% Reinforcement Areas table
reinforcementAreaTitle = Paragraph('Reinforcement Areas');
reinforcementAreaTitle.Style = {Bold(true)};
add(sec1, reinforcementAreaTitle);

reinforcementAreaTable = {...
    'Bottom steel area', s_area_bot;
    'Reinforcement ratio in tension', rho_t;
    'Top steel area', s_area_top;
    'Reinforcement ratio in compression', rho_c};
t5 = FormalTable(reinforcementAreaTable);
t5.Style = {Border('solid'), RowSep('solid'), ColSep('solid')};
add(sec1, t5);

add(ch, sec1);

% Add figures and descriptions
figTitles = {'Tension Model', 'Compression Model', 'Reinforcement Model', ...
              'M vs beta','k vs beta','Moment vs Curvature','Load vs Deflection','Stress Distribution'};
sec4 = Section('Figures');

for i = 1:length(figTitles)
    figFileName = ['figure' num2str(i) '.png'];
    saveas(figure(i), figFileName); % Save figure as PNG file
    img = Image(figFileName);
    img.Height = '2in'; % Adjust the image size to fit on one page
    img.Width = '3in';
    add(sec4, img);
    p = Paragraph(figTitles{i});
    p.Style = {OuterMargin("0pt", "0pt", "5pt", "5pt")};
    add(sec4, p);
end

add(ch, sec4);

% Add the chapter to the report
add(rpt, ch);

% Close the report
close(rpt);

% Display the location of the generated PDF
%disp(['PDF report generated: ' fileName]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Write Output file                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Create_Outputfile == 1

% Create tables for stress-strain data for tension, compression, and reinforcement models
strainT = [0, epsilon_cr, epsilon_cr * beta_1, epsilon_cr * beta_2, epsilon_cr * beta_3]';
stressT = [0, epsilon_cr * E, mu_1 * epsilon_cr * E, mu_2 * epsilon_cr * E, mu_3 * epsilon_cr * E]';
tension_table = table(strainT, stressT, 'VariableNames', {'T Strains ', 'T Stresses'});

strainC = [0, omega * epsilon_cr, epsilon_cr * lambda_cu]';
stressC = [0, epsilon_cr * E * xi * omega, epsilon_cr * E * xi * omega * mu_c]';
compression_table = table(strainC, stressC, 'VariableNames', {'epsilon_c', 'sigma_c'});

strainR = [0, kappa * epsilon_cr, epsilon_cr * chi_su]';
stressR = [0, epsilon_cr * E * n * kappa, epsilon_cr * E * n * kappa * mu_s]';
reinforcement_table = table(strainR, stressR, 'VariableNames', {'epsilon_s', 'sigma_s'});

% Write stress-strain data to Excel
writetable(tension_table, output_filename, 'Sheet', 'Tension_Data');
writetable(compression_table, output_filename, 'Sheet', 'Compression_Data');
writetable(reinforcement_table, output_filename, 'Sheet', 'Reinforcement_Data');

% Create tables for Moment, Curvature, K, delta_total, and Load_4pb
moment_curvature_table = table(Moment, Curvature, k_final, 'VariableNames', {'M', 'phi', 'k'});
delta_load_table = table(delta_total, load, 'VariableNames', {'delta', 'Load'});

% Write Moment, Curvature, K, delta_total, and Load_4pb to Excel
writetable(moment_curvature_table, output_filename, 'Sheet', 'Moment_Curvature_Data');
writetable(delta_load_table, output_filename, 'Sheet', 'Deflection_Load_Data');

% Create a table for the "Strain Values" sheet
strain_values_table = table(delta_total, load,strain_bot, est_bot, est_top, strain_top, NA_from_bot, ...
    'VariableNames', {'delta', 'Load','epsilon_bot', 'es_T_bot', 'es_T_top', 'epsilon_top', 'NA_from_bot'});

% Write the strain values to the "Strain Values" sheet
writetable(strain_values_table, output_filename, 'Sheet', 'Strain_Values');

% Transpose input parameters and add subtitles for sections
output_data = {
    'Geometry', ''; 'pointBend', pointBend; 'L', L; 'b', b; 'h', h; 'S2', S2; 'Lp', Lp; 'cover', cover; 'alpha', alpha; '', '';
    'Tension Model', ''; 'fc', fc; 'E', E; 'epsilon_cr', epsilon_cr; 'mu_1', mu_1; 'mu_2', mu_2; 'mu_3', mu_3; 'beta_1', beta_1; 'beta_2', beta_2; 'beta_3', beta_3; 'eta_1', eta_1; 'eta_2', eta_2; 'eta_3', eta_3; '', '';
    'Compression Model', ''; 'xi', xi; 'ecy', ecy; 'omega', omega; 'mu_c', mu_c; 'ecu', ecu; 'lambda_cu', lambda_cu; 'eta_c', eta_c; '', '';
    'Steel Properties', ''; 'Es', Es; 'n', n; 'kappa', kappa; 'mu_s', mu_s; 'chi_su', chi_su; 'eta_s', eta_s; '', '';
    'Reinforcement Details', ''; 'topDiameter', topDiameter; 'topCount', topCount; 'botDiameter', botDiameter; 'botCount', botCount; '', '';
    'Reinforcement Areas', ''; 's_area_bot', s_area_bot; 'rho_t', rho_t; 's_area_top', s_area_top; 'rho_c', rho_c;
};

% Convert to table and write to Excel
output_table = cell2table(output_data, 'VariableNames', {'Parameter', 'Value'});
writetable(output_table, output_filename, 'Sheet', 'Input_Parameters', 'WriteVariableNames', false);

% Create a table for the "Efficiency Data" sheet
eff_values_table = table(force_Tens, force_Comp, tot_force_rebar_bot, tot_force_rebar_top, ...
    'VariableNames', {'force_Tens', 'force_Comp', 'tot_force_rebar_bot', 'tot_force_rebar_top'});

% Write the efficiency values to the "Efficiency Data" sheet
writetable(eff_values_table, output_filename, 'Sheet', 'Efficiency_Data');

% Experimental Data of Load-Deflection if available        

% Check if experimental data is available (hasFlex = 1)
if hasFlex == 1
    % Create a table with experimental load-deflection data
    exp_data_table = table(exp_deflection_data, exp_load_data, ...
        'VariableNames', {'Deflection', 'Load'});
    
    % Write the experimental data to a new sheet
    writetable(exp_data_table, output_filename, 'Sheet', 'Experimental_Data');
else
    % If no experimental data is available, create a table that indicates this
    no_data_table = table({'Experimental data not available'}, 'VariableNames', {'Message'});
    
    % Write the message to the "Experimental Data" sheet
    writetable(no_data_table, output_filename, 'Sheet', 'Experimental_Data');
end

disp('Output data written to Excel file successfully.');
end