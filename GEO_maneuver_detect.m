function [X, maneuvercount, maneuveridx]= GEO_maneuver_detect(X, tspan, n_svs, threshold, fstep, fn, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function models maneuvers and post-maneuever trajectory
% reintegration for any number of GEO orbits.
%
% Inputs
%   X - (t x 18*N) Preliminary 2D array of state vectors (x, y, z, vx, vy, vz) for N
%       satellites, the Sun, and the Moon over t time
%   tspan - (1 x t) Vector of time points
%   n_svs - Number of satellites
%   threshold - Value to determine need for maneuver
%                   4/12/22: Max separation distance (km)
%   fstep - True anomaly spacing (rad)
%   fn - Function handle for trajectory reintegration
%   options - Tolerance options for MATLAB ODE function
%
% Outputs
%   X - (t x 18*N) 2D array of state vectors (x, y, z, vx, vy, vz) for N
%       satellites, the Sun, and the Moon over t time that accounts for maneuvers and
%       reintegration
%   maneuvercount - Number of maneuvers conducted
%   maneuveridx - Time indices of maneuvers
%
% Example
%
%   [X_twobody, maneuvercount_twobody, maneuveridx_twobody] ...
%       = GEO_maneuver_detect(X_twobody, tspan, n_svs, threshold, fstep, fn, options);
%   
%   where
%       tspan = 0:T_GEO:365*T_GEO; % sec
%       fn = @ode_J2SM_cowells;
%       options = odeset('RelTol', 1.0e-13, 'InitialStep', 1.0e-12, 'AbsTol', 1.0e-10);  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare global variables
global i_GEO r_GEO v_GEO;

% Determine need for maneuver
maneuvercount = 0;
maneuveridx = 1;
for i = 1:length(tspan)
    for j = 1:n_svs
        
        % Calculate distance between adjacent SVs
        svA = X(i, (18*(j - 1) + 1):(18*(j - 1) + 3));
        if j < n_svs
            svB = X(i, (18*(j) + 1):(18*(j) + 3));
        else
            svB = X(i, 1:3);
        end
        
        % If spacing is greater than threshold, maneuver and re-integrate
        % orbits
        if norm(svA - svB) > threshold
            
            % Maneuver all satellites in constellation for redundancy
            for k = 1:n_svs
                if k == 1
                    f = atan(svA(2)/svA(1)); % (rad) True anomaly of primary SV
                    if (svA(1) < 0 && svA(2) > 0) || (svA(1) < 0 && svA(2) < 0)
                        f = f + pi; % Correct quadrant if necessary
                    end
                else
                    f = f + fstep; % (rad) True anomaly of subsequent SVs
                end
                
                % Small routine to translate variables to actual SVs
                svidx = j + k - 1;
                if svidx > n_svs
                    svidx = svidx - n_svs;
                end
                
                % Position and velocity post-maneuver
                % Maneuver is an instantaneous radial maneuver
                r0_new = [r_GEO*cos(f)*cos(i_GEO) r_GEO*sin(f)*cos(i_GEO) r_GEO*sin(i_GEO)]';
                v0_new = [-v_GEO*sin(f)*cos(i_GEO) v_GEO*cos(f)*cos(i_GEO) v_GEO*sin(i_GEO)]';
                
                % New initial state
                X0_new = [r0_new; v0_new; X(i, (18*(svidx - 1) + 7):(18*(svidx - 1) + 18))'];
                
                % New time span
                tspan_new = tspan(i:end);
                
                % Integrate post-maneuver orbit
                [~, X(i:end, (18*(svidx - 1) + 1):(18*(svidx - 1) + 18))] ...
                    = ode45(fn, tspan_new, X0_new, options);
            end
            
            % Keep record of maneuvers
            maneuvercount = maneuvercount + 1;
            maneuveridx(end + 1) = i;
            continue;
        end
    end
end
maneuveridx(end + 1) = length(tspan); % Add end of time span
end

            
%% BACKUP

% % Calculate distance between SV A and SV B
% distOverTime = zeros(length(tspan), n_svs);
% for i = 1:length(tspan)
%     for j = 1:n_svs
%         svA = X(i, (6*(j - 1) + 1):(6*(j - 1) + 3));
%         if j < n_svs
%             svB = X(i, (6*(j) + 1):(6*(j) + 3));
%         else
%             svB = X(i, 1:3);
%         end
%         distOverTime(i, j) = norm(svA - svB);
%     end
% end

% Calculate orbital radius over time
% rOverTime = zeros(length(tspan), n_svs);
% for i = 1:length(tspan)
%     for j = 1:n_svs
%         rOverTime(i, j) = norm(X(i, (6*(j - 1) + 1):(6*(j - 1) + 3)));
%     end
% end

% for j = 1:n_svs
% %     maneuveridx_list{j}(end + 1) = 1;
%     maneuveridx = find((abs(r_GEO - rOverTime(:, j)) > threshold) == 1, 1);
%     while ~isempty(maneuveridx)
%         
%         % Keep record of maneuvers
%         maneuveridx_list{j}(end + 1) = maneuveridx;
%         
%         % Execute maneuver
%         r0_new = X(maneuveridx, (6*(j - 1) + 1):(6*(j - 1) + 3))';
%         r0_new(2) = sqrt(r_GEO^2 - X(maneuveridx, (6*(j - 1) + 1))^2); % Change y-position only
%         if X(maneuveridx-1, (6*(j - 1) + 2)) < 0
%             r0_new(2) = -r0_new(2); % Adjust sign if necessary
%         end
%         f_new = atand(r0_new(2)/r0_new(1)); %Calculate new true anomaly
%         v0_new = [-v_GEO*sind(f_new)*cos(i_GEO) v_GEO*cosd(f_new)*cosd(i_GEO) v_GEO*sind(i_GEO)]';
%         
%         % New initial state
%         X0_new = [r0_new; v0_new];
%         
%         % New time span
%         tspan_new = tspan(maneuveridx:end);
%         
%         % Integrate 3D two-body problem
%         [t, X(maneuveridx:end, (6*(j - 1) + 1):(6*(j - 1) + 6))] = ode45(fn, tspan_new, X0_new, options);
%         
%         % Calculate new orbital radii over time
%         rOverTime(maneuveridx:end, j) = zeros(length(tspan_new), 1);
%         for k = maneuveridx:length(tspan)
%             rOverTime(k, j) = norm(X(k, (6*(j - 1) + 1):(6*(j - 1) + 3)));
%         end
%         
%         maneuvercount(j) = maneuvercount(j) + 1; % Iterate count
%         maneuveridx = find((abs(r_GEO - rOverTime(:, j)) > threshold) == 1, 1); % Find next need for maneuver
%         
%     end
%     maneuveridx_list{j}(end + 1) = length(tspan); % Add end of time span
% end
% 
% 
% 
% 
% % % Determine need for maneuver per SV
% % maneuvercount = zeros(1, n_svs);
% % manueveridx_list = cell(1, n_svs);
% % for k = 1:n_svs
% %     maneuveridx_list{k}(end + 1) = 1;
% %     maneuveridx = find((distOverTime(:, k) > threshold) == 1, 1);
% %     while ~isempty(maneuveridx)
% %         
% %         % Keep record of maneuvers
% %         maneuveridx_list{k}(end + 1) = maneuveridx;
% %         
% %         % Execute maneuver
% %         r = X(maneuveridx, (6*(k - 1) + 1):(6*(k - 1) + 3))';
% %         f = atand(r(2), r(1));
% %         r0_new = [r_GEO*cosd(f) r_GEO*sind(f) r(3)]'; % maneuver radially instantaneously
% % %         r0_new(2) = r_GEO*sind(f); % maneuver radially instantaneously
% % %         r0_new(2) = sqrt(r_GEO^2 - X(maneuveridx, 6*(k - 1) + 1)^2); % Change y-position only
% % %         if X(maneuveridx-1, 6*(k - 1) + 2) < 0
% % %             r0_new(2) = -r0_new(2); % Adjust sign if necessary
% % %         end
% % %         f_new = atand(r0_new(2)/r0_new(1)); %Calculate new true anomaly
% %         v0_new = [-v_GEO*sind(f)*cos(i_GEO) v_GEO*cosd(f)*cosd(i_GEO) v_GEO*sind(i_GEO)]';
% %         
% %         % TODO: Check if adjacent SV needs to maneuver as well to maintain
% %         % separation distance
% %         svA = r0_new;
% %         svB = X(maneuveridx, (6*(k) + 1):(6*(k) + 3))';
% %         if norm(svA - svB) > threshold
% %             svB = [r_GEO*cosd(f+fstep) r_GEO*sind(f+step) svB(3)]; % maneuver radially instantaneously
% %         end
% %         
% %         % New initial state
% %         X0_new = [r0_new; v0_new];
% %         
% %         % New time span
% %         tspan_new = tspan(maneuveridx:end);
% %         
% %         % Integrate 3D two-body problem
% %         [t, X(maneuveridx:end, (6*(k - 1) + 1):(6*(k - 1) + 6))] ...
% %             = ode45(@fn, tspan_new, X0_new, options);
% %         
% %         % Calculate new orbital radii over time
% % %         rOverTime(maneuveridx:end) = zeros(length(tspan_new), 1);
% % %         for k = maneuveridx:length(tspan)
% % %             rOverTime(k) = norm(X_twobody(k, 1:3));
% % %         end
% %         
% %         % TODO: Calculate distance between SV A and SV B
% %         distOverTime(maneuveridx:end, k) = zeros(length(tspan), n_svs);
% %         for i = 1:length(tspan)
% %             for j = 1:n_svs
% %                 svA = X(i, (6*(j - 1) + 1):(6*(j - 1) + 3));
% %                 if j < n_svs
% %                     svB = X(i, (6*(j) + 1):(6*(j) + 3));
% %                 else
% %                     svB = X(i, 1:3);
% %                 end
% %                 distOverTime(i, j) = norm(svA - svB);
% %             end
% %         end
% %         
% %         maneuvercount = maneuvercount + 1; % Iterate count
% %         maneuveridx = find((rOverTime < threshold) == 1, 1); % Find next need for maneuver
% %         
% %     end
% %     maneuveridx_list(end + 1) = length(tspan); % Add end of time span
% % end

