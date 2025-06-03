function [Envelope,T,C,R,RC] = calculateEnvelope_new_2(kappa, omega, epsilon_cr,beta_all, k111, M111, ...
                                      k211, M211, k212, M212, k221, M221, k222, M222, ...
                                      k311, M311, k312, M312, k321, M321, k322, M322, ...
                                      k411, M411, k412, M412, k421, M421, k422, M422, ...
                                      k4222,M4222,beta_1,beta_2,beta_3,alpha,T,C,R,RC)

    % persistent warningDisplayed; % Persistent variable to track if the warning has been displayed
    % 
    % if isempty(warningDisplayed) % Check if the warning has been displayed before
    % warningDisplayed = false;
    % end
% CALCULATEENVELOPE_NEW
% Example skeleton for how you might select stage-based k, M arrays
% and update T (tension zone), C (compression zone), R (steel zone).
%
% INPUTS (abbreviated):
%   kappa, omega, epsilon_cr : strain/curvature parameters
%   beta_all                 : N×1 vector, e.g. your incremental parameter
%   beta_z1..z4              : thresholds for tension zones
%   k111, M111, etc.         : each is an N×1 vector for that stage
%   alpha                    : geometric offset for strain calculations
%
% OUTPUT:
%   Envelope : N×2 array, [k(i), M(i)] for each increment.
%
% NOTES:
% 1) T, C, R are integers indicating the zone: T = 1..4, C = 1..2, R = 1..2
% 2) You have 13 or 16 “feasible” stages. Adjust the switch statement as needed.
% 3) The code below shows how to pick k/M based on [T,C,R], then update T,C,R.

% Initialize
n = size(beta_all, 1);
Envelope = zeros(n,2);   % store [k(i), M(i)] at each step
econ_y   = zeros(n,1);  % store compression fiber strain
est_y    = zeros(n,1);  % store steel strain
esc_y    = zeros(n,1);  % store steel strain


% Loop over each increment of beta_all
for i = 1 : n
    %---------------------------------------------
    % 1) Determine which stage:  e.g. '111', '212', etc.
    %---------------------------------------------
    stageStr = [num2str(T), num2str(C), num2str(R), num2str(RC)];

    %---------------------------------------------
    % 2) Select the correct k, M arrays for that stage
    %---------------------------------------------
    switch stageStr
        case '1111'
            kStage = k111; MStage = M111;
        case '2111'
            kStage = k211; MStage = M211;
        case '2121'
            kStage = k212; MStage = M212;
        case '2211'
            kStage = k221; MStage = M221;
        case '2221'
            kStage = k222; MStage = M222;
        case '3111'
            kStage = k311; MStage = M311;
        case '3121'
            kStage = k312; MStage = M312;
        case '3211'
            kStage = k321; MStage = M321;
        case '3221'
            kStage = k322; MStage = M322;
        case '4111'
            kStage = k411; MStage = M411;
        case '4121'
            kStage = k412; MStage = M412;
        case '4211'
            kStage = k421; MStage = M421;
        case '4221'
            kStage = k422; MStage = M422;
        case '4222'
            kStage = k4222; MStage = M4222;    
        otherwise
            warning('Stage "%s" not recognized', stageStr);
    end

    %---------------------------------------------
    % 3) Now get the "current" k(i) and M(i) for this increment
    %---------------------------------------------
    Envelope(i,1) = kStage(i,1);  % k(i)
    Envelope(i,2) = MStage(i,1);  % M(i)

    %---------------------------------------------
    % 4) Strain calculations at top and steel, etc.
    %    (using your formula for econ_y, est_y)
    %---------------------------------------------
    k_i = kStage(i);
    beta_i = beta_all(i);

    % econ_y = top fiber strain
    econ_y(i) = k_i * beta_i * epsilon_cr / (1 - k_i);

    % est_y = steel tension strain
    est_y(i)  = ((-alpha + k_i) * beta_i * epsilon_cr) / (k_i - 1);
    
    % esc_y = compression steel strain
    esc_y(i)  = abs(((k_i-1+alpha) * beta_i * epsilon_cr) / (k_i - 1));

    %---------------------------------------------
    % 5) Update zones T, C, R by checking thresholds
    %---------------------------------------------
    % T: tension zone in concrete could go from 1->2->3->4
    if (T == 1) && (beta_i >= 1)
        T = 2;
    elseif (T == 2) && (beta_i >= beta_1)
        T = 3;
    elseif (T == 3) && (beta_i >= beta_2)
        T = 4;
    end

    % C: compression zone in concrete (only 1->2)
    if (C == 1) && (econ_y(i) >= (omega * epsilon_cr))
        C = 2;
    end

    % R: rebar zone (steel) 1->2 if yield
    if (R == 1) && (est_y(i) >= kappa * epsilon_cr)
        R = 2;
    end

    % RC: compression steel yielding check
    if (RC == 1) && (esc_y(i) >= kappa * epsilon_cr) && strcmp(stageStr, '4221') || strcmp(stageStr,'4222')
        RC = 2;
    end

    % if (RC == 1) && (esc_y(i) >= kappa * epsilon_cr) && ~warningDisplayed
    %     warning('Compression Steel assumed Elastic because it yeilds while tension steel or conctete in compression is still elastic.Reduce Tension Steel');
    % warningDisplayed = true; % Set the flag to true to prevent further warnings
    % end

end % end for

end % function
