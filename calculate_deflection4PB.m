function delta_total = calculate_deflection4PB(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, Lp2, Mmax,S2)
    % Initialize variables
    delta_total = zeros(length(mom), 1);
    first_delta_total = []; % To store the delta_total value at the first i meeting the condition
    last_delta_total = 0;


    for i = 1:length(mom)
        if mom(i) <= M_cr && cv(i) < Cmax && cv(i) <= Phi_cr
            %delta_total(i, 1) = 0.5 * cv(i) * (L/2 - S2/2) * (L/3 - S2/3) + cv(i) * S2 * (L/2 - S2/4) / 2;
            delta_total(i, 1) = (0.5 * cv(i) * (L/2 - S2/2) * (2/3*(L/2 - S2/2))) + (cv(i) * S2/2 * (L/2 - S2/4));
            last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value

        elseif mom(i) < M_cr && cv(i) <= Cmax && cv(i) > Phi_cr
            x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp/2));
            aa = x_cv_cr;
            bb = (L/2)-(Lp/2);
            cc = L/2;
            if aa>=bb
                aa=bb;
            end


            x1 = (2/3) * aa;
            x2 = ((bb - aa) / 2) + aa;
            x3 = aa + ((2/3) * (bb - aa));
            x4 = bb + ((cc - bb) / 2);

            A1 = 0.5 * aa * Phi_cr;
            A2 = Phi_cr * (bb - aa);
            A3 = 0.5 * (bb - aa) * (cv(i) - Phi_cr);
            A4 = (cc - bb) * cv(i);

            delta1 = A1 * x1;
            delta2 = A2 * x2;
            % if delta2 <= 0
            %     delta2 = 0;
            % end
            delta3 = A3 * x3;
            %  if delta3 <= 0
            %     delta3 = 0;
            % end
            delta4 = A4 * x4;

            delta_total(i, 1) = delta1 + delta2 + delta3 + delta4;
            last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value

        elseif mom(i) > M_cr && cv(i) <= Cmax
            x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp/2));
            aa = x_cv_cr;
            bb = (L/2)-(Lp/2);
            cc = L/2;
            if aa>=bb
               aa=bb;
            end


            x1 = (2/3) * aa;
            x2 = ((bb - aa) / 2) + aa;
            x3 = aa + ((2/3) * (bb - aa));
            x4 = bb + ((cc - bb) / 2);

            A1 = 0.5 * aa * Phi_cr;
            A2 = Phi_cr * (bb - aa);
            A3 = 0.5 * (bb - aa) * (cv(i) - Phi_cr);
            A4 = (cc - bb) * cv(i);

            delta1 = A1 * x1;
            delta2 = A2 * x2;
            delta3 = A3 * x3;
            delta4 = A4 * x4;

            delta_total(i, 1) = delta1 + delta2 + delta3 + delta4;
            last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value

        elseif mom(i) < Mmax && cv(i) > Cmax
            aa = L/2 - (Lp2/2);
            bb = L/2 - aa;
            x1 = (2/3) * aa;
            x2 = aa + (bb / 2);

            if aa>=bb
               aa=bb;
            end


            delta1 = x1 * (0.5 * Phi_cr * aa);
            delta2 = x2 * bb * cv(i);
            delta_total(i, 1) = delta1 + delta2;

            % Record the first delta_total value
            if isempty(first_delta_total)
                first_delta_total = delta_total(i, 1);
            end
            delta_total(i, 1) = (last_delta_total - first_delta_total) + delta1 + delta2;
        end
    end
end






% function delta_total = calculate_deflection(mom, cv, M_cr, Cmax, Phi_cr, L, Lp, Lp2, Mmax,S2)
%     % Initialize variables
%     delta_total = zeros(length(mom), 1);
%     first_delta_total = []; % To store the delta_total value at the first i meeting the condition
%     last_delta_total = 0;
% 
% 
%         if Cmax <= Phi_cr
%             result = cv(cv >= Phi_cr);
%             num_points = length(result);  % Find the number of points in the result vector
%             Lp_fact = linspace(1, 1, num_points);  % Create a linspace vector of the same size
%             Lp_vect = Lp_fact*Lp;
%         else
%             result0 = cv(cv <= Phi_cr);
%             result = cv(cv >= Phi_cr & cv <= Cmax);
%             result1 = cv(cv >= Cmax);
% 
%             num_points0 = length(result0);
%             num_points = length(result);
%             num_points1 = length(result1);  % Find the number of points in the result vector
%             Lp_fact0 = linspace(0, 0, num_points0);  % Create a linspace vector of the same size
%             Lp_fact = linspace(0, 1, num_points);  % Create a linspace vector of the same size
%             Lp_fact1 = linspace(1, 1, num_points1);  % Create a linspace vector of the same size
%             Lp_vect = [(Lp_fact0*Lp),(Lp_fact*Lp),(Lp_fact1*Lp)];
%         end
% 
%     for i = 1:length(mom)
%         if mom(i) <= M_cr && cv(i) < Cmax && cv(i) <= Phi_cr
%             delta_total(i, 1) = 0.5 * cv(i) * (L/2 - S2/2) * (L/3 - S2/3) + cv(i) * S2 * (L/2 - S2/4) / 2;
%             last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value
% 
%         elseif mom(i) < M_cr && cv(i) <= Cmax && cv(i) > Phi_cr
% 
%                 if Lp_vect(i)>= S2
%                 x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp_vect(i)/2));
%                 else
%                 x_cv_cr = (M_cr / mom(i)) * ((L/2)-(S2/2));
%                 end
% 
%             %x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp/2));
%             aa = x_cv_cr;
%             bb = (L/2)-(Lp_vect(i)/2);
%             cc = L/2;
% 
%             x1 = (2/3) * aa;
%             x2 = ((bb - aa) / 2) + aa;
%             x3 = aa + ((2/3) * (bb - aa));
%             x4 = bb + ((cc - bb) / 2);
% 
%             A1 = 0.5 * aa * Phi_cr;
%             A2 = Phi_cr * (bb - aa);
%             A3 = 0.5 * (bb - aa) * (cv(i) - Phi_cr);
%             A4 = (cc - bb) * cv(i);
% 
%             delta1 = A1 * x1;
%             delta2 = A2 * x2;
%             delta3 = A3 * x3;
%             delta4 = A4 * x4;
% 
%             delta_total(i, 1) = delta1 + delta2 + delta3 + delta4;
%             last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value
% 
%         elseif mom(i) > M_cr && cv(i) <= Cmax
%            % x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp/2));
% 
% 
%                 if Lp_vect(i)>= S2
%                 x_cv_cr = (M_cr / mom(i)) * ((L/2)-(Lp_vect(i)/2));
%                 else
%                 x_cv_cr = (M_cr / mom(i)) * ((L/2)-(S2/2));
%                 end
% 
%             aa = x_cv_cr;
%             bb = (L/2)-(Lp_vect(i)/2);
%             cc = L/2;
% 
%             x1 = (2/3) * aa;
%             x2 = ((bb - aa) / 2) + aa;
%             x3 = aa + ((2/3) * (bb - aa));
%             x4 = bb + ((cc - bb) / 2);
% 
%             A1 = 0.5 * aa * Phi_cr;
%             A2 = Phi_cr * (bb - aa);
%             A3 = 0.5 * (bb - aa) * (cv(i) - Phi_cr);
%             A4 = (cc - bb) * cv(i);
% 
%             delta1 = A1 * x1;
%             delta2 = A2 * x2;
%             delta3 = A3 * x3;
%             delta4 = A4 * x4;
% 
%             delta_total(i, 1) = delta1 + delta2 + delta3 + delta4;
%             last_delta_total = delta_total(i, 1); % Update the last_delta_total with the current value
% 
%         elseif mom(i) < Mmax && cv(i) > Cmax
%             aa = L/2 - (Lp2/2);
%             bb = L/2 - aa;
%             x1 = (2/3) * aa;
%             x2 = aa + (bb / 2);
% 
%             delta1 = x1 * (0.5 * Phi_cr * aa);
%             delta2 = x2 * bb * cv(i);
%             delta_total(i, 1) = delta1 + delta2;
% 
%             % Record the first delta_total value
%             if isempty(first_delta_total)
%                 first_delta_total = delta_total(i, 1);
%             end
%             delta_total(i, 1) = (last_delta_total - first_delta_total) + delta1 + delta2;
%         end
%     end
%end