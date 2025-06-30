classdef Control_SMC_MPC
    properties(SetAccess = private)
        % AUV parameters
        dT %--> Controller sampling time (1x1 double, unit: sec.)
        Tlim %--> Thruster output limit (1x1 double, unit: N)

        % Dynamic properties
        taw %--> Control thrust (1x8 array of double, unit: N)
        n %--> Global pose feedback (1x6 array of double, unit: m and/or rad)
        v %--> Local velocity estimate (1x6 array of double, unit: m/s and/or rad/s)

        % HSMC
        eta_tilde_prev % Previous tracking error
        t % Time tracking for control
        Sd % Initial sliding surface (constant)
        S_eta_int % Integral of sliding error
        tb
        j0
        j1
        s_prev

        % % AUV parameters
        % m = 25.3050;      % Mass (kg)
        % Iz = 0.0;      % Yaw inertia (kg*m^2)
        % Xu_dot = 360.3;
        % Yv_dot = 360.2;
        % Nr_dot = 8.55;
        % Xu = 0;     % Surge drag coefficient
        % Yv = 0;     % Sway drag coefficient
        % Nr = 0;      % Yaw drag coefficient
        % Du = 298.2797;
        % Dv = 299.8409;
        % Dr = 2.8145;

        % BLUE AUV parameters
        m = 13.5;      % Mass (kg)
        Iz = 0.37;      % Yaw inertia (kg*m^2)
        Xu_dot = 6.357;
        Yv_dot = 7.121;
        Nr_dot = 0.2215;
        Xu = 13.7;     % Surge drag coefficient
        Yv = 0;     % Sway drag coefficient
        Nr = 0;      % Yaw drag coefficient
        Du = 141;
        Dv = 217;
        Dr = 1.5;
        
        % MPC parameters
        T = 1.5;        % Prediction horizon (s)
        N;         % Control intervals
        
        % Input constraints
        Tx_max = 50.0;  % Max thrust (N)
        Ty_max = 5.0;
        Mz_max = 50.0;  % Max moment (Nm)
        args = struct;
        
        % Weighting matrices
        Q = diag([50.0, 50.0, 5.0, 0.1, 0.1, 0.1]);  % State weights
        R = diag([0.01, 0.01, 0.01]);                % Input weights
        P

        % Backstepping gain (for path tracking)
        k_r = 1.0;       % Heading convergence rate
        prev_delta = 0;
        prev_beta = 0;
        path_spline;
        max_s=1;
        
        % Path tracking parameters
        path_points = [];
        path_Z =[];
        idx = 1;
        path_plot
        auv_marker;
        lookahead = 0.8;  % set your desired threshold value here
        current_s = 0;
        
        % CasADi objects
        solver;         % NLP solver
        F_rk4;          % RK4 integrator
        figure_handle
    end
    properties
        alpha_0  
        alpha_c
        delta 
        K_d 
        K_i 
        seg
        vel
    end
    methods
        function obj = Control_SMC_MPC(sens, dT)
            % Constructor
            obj.taw = zeros(1,8);
            obj.n = [sens.orbX, sens.orbY, sens.ps, sens.imu_euler(1), sens.imu_euler(2), sens.imu_euler(3)];
            obj.v = zeros(1,6);
            obj.dT = dT;
            obj.Tlim = 10;

            % HSMC
            % Initialize Controller States
            obj.t = 0;
            obj.S_eta_int = zeros(3,1);
            obj.Sd = zeros(3,1); % This will be initialized when first computed
            obj.j0 = [];
            obj.j1 = [];
            obj.alpha_0 = 500;  
            obj.alpha_c =1;
            obj.delta = 0.1;
            obj.K_d = diag([100, 5, 5]); % Gain for error correction
            obj.K_i = diag([0.001, 0.0001, 0.0001]); % Integral gain
            obj.vel=0.5;
            obj.seg=1;
            obj.eta_tilde_prev = [0;0;0];
            
            % MPC
            obj.N = ceil(obj.T / dT); % Adjust control intervals
            obj.s_prev = [0,0,0];
            % Constructor: Set up CasADi and MPC solver.
            obj = obj.setup_casadi();   
        end

        function obj = setup_casadi(obj)
            % Initialize CasADi symbolic variables and solver.
            import casadi.*
        
            % States (x, y, psi, u, v, r) and inputs (Tx, Ty, Mz)
            x = SX.sym('x', 6);
            u = SX.sym('u', 3);
        
            n_states = length(x);
            n_controls = length(u);
        
            % Extract states
            psi = x(3);
            v_body = x(4:6);  % [u; v; r]
        
            % Rotation matrix R(psi)
            Rot = [cos(psi), -sin(psi);
                   sin(psi),  cos(psi)];
        
            % Kinematics: [ẋ, ẏ, ψ̇]
            eta_dot = Rot * v_body(1:2);
            psi_dot = v_body(3);
        
            % Inertia matrix (diagonal)
            Mu = obj.m + obj.Xu_dot;
            Mv = obj.m + obj.Yv_dot;
            Mr = obj.Iz + obj.Nr_dot;
            M = diag([Mu, Mv, Mr]);
        
            % Coriolis-centripetal matrix C(v)
            C = [  0,   0, -Mv * v_body(2);
                   0,   0,  Mu * v_body(1);
                 Mv * v_body(2), -Mu * v_body(1), 0];
        
            % Damping matrix D(v) = linear + quadratic drag
            D = diag([-obj.Xu * v_body(1) + obj.Du * abs(v_body(1));
                      -obj.Yv * v_body(2) + obj.Dv * abs(v_body(2));
                      -obj.Nr * v_body(3) + obj.Dr * abs(v_body(3))]);
        
            % Dynamics: [u̇, v̇, ṙ]
            v_dot = M \ (u - C * v_body - D * v_body);
        
            % Full state derivative
            xdot = [eta_dot; psi_dot; v_dot];
        
            % RK4 integrator (discrete-time)
            k1 = xdot;
            k2 = xdot + obj.dT/2 * k1;
            k3 = xdot + obj.dT/2 * k2;
            k4 = xdot + obj.dT * k3;
            x_next = x + obj.dT/6 * (k1 + 2*k2 + 2*k3 + k4);
            obj.F_rk4 = Function('F_rk4', {x, u}, {x_next});
        
            % Multiple shooting: Decision variables
            X_mpc = SX.sym('X_mpc', 6, obj.N + 1);  % State trajectory
            U_mpc = SX.sym('U_mpc', 3, obj.N);      % Control inputs
        
            % Parameters: Initial state and reference
            x_init = SX.sym('x_init', 6);
            x_ref = SX.sym('x_ref', 6);
        
            %% Cost and constraints
            cost = 0;
            g = [];
        
            % Initial state constraint
            g = [g; X_mpc(:, 1) - x_init];
        
            % Loop over horizon
            for i = 1:obj.N
                % State error cost
                error = X_mpc(:, i) - x_ref;
                cost = cost + error' * obj.Q * error;
        
                % Control cost
                cost = cost + U_mpc(:, i)' * obj.R * U_mpc(:, i);
        
                % Continuity constraint (multiple shooting)
                x_next = obj.F_rk4(X_mpc(:, i), U_mpc(:, i));
                g = [g; X_mpc(:, i + 1) - x_next];
            end
        
            % === Terminal cost using Lyapunov function ===
            % You can tune P or compute it analytically if needed
            obj.P = 10 * obj.Q;  % Make sure this is set before calling this function
            x_terminal = X_mpc(:, end);
            terminal_error = x_terminal - x_ref;
            terminal_cost = terminal_error' * obj.P * terminal_error;
            cost = cost + terminal_cost;
        
            % === Bounds for constraints ===
            obj.args.lbg = zeros(n_states*(obj.N+1),1);
            obj.args.ubg = zeros(n_states*(obj.N+1),1);
        
            % Define bounds on U_mpc and X_mpc separately
            lbx_u = repmat([-obj.Tx_max; -obj.Ty_max; -obj.Mz_max], obj.N, 1);
            ubx_u = repmat([ obj.Tx_max;  obj.Ty_max;  obj.Mz_max], obj.N, 1);
        
            % For X_mpc (no bounds here, so just use -inf/inf)
            lbx_x = -inf * ones(6 * (obj.N + 1), 1);
            ubx_x =  inf * ones(6 * (obj.N + 1), 1);
        
            % Concatenate
            obj.args.lbx = [lbx_u; lbx_x];
            obj.args.ubx = [ubx_u; ubx_x];
        
            % NLP problem
            OPT_variables = [reshape(U_mpc,n_controls*obj.N,1);reshape(X_mpc,n_states*(obj.N+1),1)];
            nlp = struct(...
                'f', cost, ...                  % Cost function
                'x', OPT_variables, ...        % Decision variables
                'g', g, ...                     % Constraints
                'p', [x_init; x_ref] ...       % Parameters
            );
        
            % Solver options
            opts = struct;
            opts.ipopt.max_iter = 2000;
            opts.ipopt.print_level = 0;
            opts.print_time = 0;
            opts.ipopt.acceptable_tol = 1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
        
            % Create solver
            obj.solver = nlpsol('solver', 'ipopt', nlp, opts);
        end


        function [u_opt, X_pred] = solve_mpc(obj, x0, x_ref)
            import casadi.*
            
            % Initial guess
            X_guess = repmat(x0, 1, obj.N+1);
            U_guess = zeros(3, obj.N);
            
            % Solve NLP
            sol = obj.solver(...
                'x0', [U_guess(:); X_guess(:)], ...
                'lbg', obj.args.lbg, ...
                'ubg', obj.args.ubg, ...
                'lbx', obj.args.lbx, ...
                'ubx', obj.args.ubx, ...
                'p', [x0(:); x_ref(:)]);
            
            % Extract solution
            U_opt = reshape(full(sol.x(1:3*obj.N)), 3, obj.N);
            X_pred = reshape(full(sol.x(3*obj.N+1:end)), 6, obj.N+1);
            u_opt = U_opt(:,1);
        end

        function [obj, f_ver] = Stabilize(obj, eta_real,eta_desired)
            % Apply Model-Free HOSMC for Stabilization in Z, Roll, Pitch
            eta_dot_d=[0;0;0];

            % Compute control thrust using HOSMC
            [obj, u] = obj.model_free_control(eta_real, eta_desired,eta_dot_d);
      
            % Distribute thrust to vertical thrusters
            j= [ 1, 1, 1, 1;
                 0.235, -0.235, -0.235, 0.235;
                 -0.1400, -0.1400, 0.1400, 0.1400;];

            j_inv = pinv(j);

            % Apply thrust limits
            f_ver = j_inv * u;
            f_ver = min(max(f_ver, -obj.Tlim), obj.Tlim);
        end

        function [obj, u] = model_free_control(obj, eta, eta_d,eta_dot_d)
            % Model-Free High-Order Sliding Mode Controller for Z, Roll, Pitch
            
            t0 = 0; % Initial time

            % Compute Tracking Error
            eta_tilde = eta - eta_d;
            [obj,eta_dot] = obj.exact_differentiator(eta_tilde,obj.eta_tilde_prev);
            obj.eta_tilde_prev = eta_tilde;
            eta_dot_tilde = eta_dot - eta_dot_d; % Estimated velocity error

            if obj.t == 0.0
                obj.tb = abs(eta_tilde(1))/obj.vel;
            end

            % Compute Time-Base Generator
            % if obj.t <= obj.tb
            %     denominator = max(obj.tb - t0, 1e-6); % Prevent division by zero
            %     xi = 10 * ((obj.t - t0) / denominator)^3 - ...
            %          15 * ((obj.t - t0) / denominator)^4 + ...
            %           6 * ((obj.t - t0) / denominator)^5;
            %     xi_dot = 30 * ((obj.t - t0)^2 / denominator^3) - ...
            %              60 * ((obj.t - t0)^3 / denominator^4) + ...
            %              30 * ((obj.t - t0)^4 / denominator^5);
            % 
            %     alpha = obj.alpha_0 * (xi_dot / ((1 - xi) + obj.delta));
            % else
            %     alpha = obj.alpha_c;
            % end
 
            % Compute time-varying gain alpha
            if obj.t <= obj.tb

                denominator = max(obj.seg * obj.tb - t0, 1e-6); % Prevent division by zero

                alpha = obj.alpha_c + (obj.alpha_0 - obj.alpha_c) * (1 - ((obj.t - t0) / denominator)^2);
            else
                alpha = obj.alpha_c;
            end

            % Compute Sliding Surface
            S = eta_dot_tilde + alpha .* eta_tilde;

            % Keep `Sd` Constant (Initialize Only Once)
            if all(obj.t == 0.0)
                obj.Sd = S; % Store initial value and never change it
            end
            
            % Compute Sliding Error & Extended Error
            S_d = obj.Sd / (1 +1000 * obj.t);
            % S_d = obj.Sd * exp(-2*obj.t); % Exponential decay with constant Sd
            S_eta = S - S_d;
            obj.S_eta_int = obj.S_eta_int + sign(S_eta) .* obj.dT; % Integral update
            S_r = S_eta + obj.K_i * obj.S_eta_int; % Extended error

            % Compute Control Law
            tau_eta = - obj.K_d * S_r;
            u = tau_eta;

            % Update Time
            obj.t = obj.t + obj.dT;
        end

        function [obj,derivative] = exact_differentiator(obj,reading,prev_reading)
            % Exact Differentiator (Super-Twisting Algorithm)
            % Inputs:
            %   reading    - Current signal value
            %   dt         - Sampling time
            % Output:
            %   derivative - Estimated derivative of the signal
            
            % Parameters (Tuning gains for sliding mode)
            % lambda1 = 0.5;  
            % lambda2 = 0.1; 
            % 
            % if isempty(obj.j0)
            %     % Initialize on first function call
            %     obj.j0 = 0;  
            %     obj.j1 = 0;        
            % end
            % 
            % % Compute error
            % e = reading - obj.j0;
            % 
            % % Super-Twisting update equations
            % obj.j0 = obj.j0 + (obj.j1 + sqrt(abs(e)) .* sign(e) .* lambda1 ) .* obj.dT;
            % obj.j1 = obj.j1 + (sign(e).*lambda2) .* obj.dT;
            % 
            % % Output estimated derivative
            % derivative = obj.j1;
            derivative = (reading-prev_reading)/obj.dT;
        end

        function [obj, f_hor] = Track(obj, pose_des, pose_curr)
            %{ 
              Return tracking thrusts for 4 horizontal thrusters
              """
              Inputs:
                    obj: class instance
                        class object
                    pose_des: 1*3 double array
                        desired pose (x, y and yaw)
                    pose_curr : 1*3 double
                        current pose (x, y and yaw)
              """
              Outputs:
                    obj: class instance
                        class object
                    f_hor: 1x4 double array
                        horizontal thrust values
               """
            %}

            % Initial state and reference (current state and desired state)
            [u, ~] = obj.solve_mpc(pose_curr, pose_des);

            % Distribute thrust to vertical thrusters
            Tht = pi/4; %--> Thruster mounting angle
            j=[ cos(Tht) , cos(Tht), cos(Tht) , cos(Tht);
                -sin(Tht), sin(Tht), -sin(Tht), sin(Tht);
                -0.375   , 0.375   , 0.375    , -0.375];
         
            j_inv = pinv(j);

            % Apply thrust limits
            f_hor = j_inv * u;
            f_hor = min(max(f_hor, -obj.Tlim), obj.Tlim);     
        end

        function obj = set_path(obj, path3D)
            % Create parametric spline
            % Separate XY and Z components
            obj.path_points = path3D(:, 1:2);  % X, Y
            obj.path_Z = path3D(:, 3);         % Z
            K = linspace(0, 1, size(obj.path_points,1))';
            obj.path_spline = spline(K, obj.path_points');
            obj.max_s = K(end);
            
            % Initialize 3D visualization
            obj.figure_handle = figure;
            figure(obj.figure_handle);
            clf;
            obj.path_plot = plot3(obj.path_points(:,1), obj.path_points(:,2), obj.path_Z, 'b-', 'LineWidth', 1.5);
            hold on;
            % After plotting the path
            obj.auv_marker = plot3(0, 0, 0, 'ro', 'MarkerSize', 8, ...
                'MarkerFaceColor', 'r', 'DisplayName', 'AUV');
            hold on;
            legend('Path', 'AUV');
            xlabel('East (m)');
            ylabel('North (m)');
            zlabel('Depth (m)');
            title('AUV Path Tracking (3D)');
            grid on;
            axis equal;
            view(3);  % 3D view
        end
        
        function [obj, f_hor] = track_path(obj, x_current)
            % 1. Find lookahead point on path
            [~, s] = obj.find_closest_point(x_current);
            [lookahead_pt, ~] = obj.find_lookahead_point(x_current);
            [~, obj.idx] = min(vecnorm(obj.path_points - lookahead_pt(1:2)', 2, 2));
            obj.current_s = s; % Update path parameter
            
            % 2. Compute path properties at lookahead point
            tangent = ppval(fnder(obj.path_spline,1), s);
            tangent = tangent/(sqrt((tangent(1))^2+(tangent(2))^2));
            psi_path = atan2(tangent(2), tangent(1));
            
            % 3. Calculate cross-track error (y1) and heading error (psi_e)
            normal = [-tangent(2); tangent(1)]; % Path normal vector
            rel_pos = x_current(1:2)' - lookahead_pt(:);
            y1 = normal' * rel_pos;
            psi_e = wrapToPi(x_current(3) - psi_path);
            
            % 4. Compute path curvature at lookahead point
            ddpt = ppval(fnder(obj.path_spline,2), s);
            c_c = (tangent(1)*ddpt(2) - tangent(2)*ddpt(1)) / (sqrt((tangent(1))^2+(tangent(2))^2))^3;
            
            % 5. Line-of-Sight (LOS) guidance
            delt = -atan2(y1, obj.lookahead); % LOS angle
            
            % 6. Backstepping control law (Eq.13 from paper)
            beta = atan2(x_current(5), x_current(4)); % Sideslip angle
            
            % Numerical derivatives
            if isempty(obj.prev_beta)
                obj.prev_beta = beta;
                obj.prev_delta = delt;
            end
            beta_dot = (beta - obj.prev_beta)/obj.dT;
            delta_dot = (delt - obj.prev_delta)/obj.dT;
            obj.prev_beta = beta;
            obj.prev_delta = delt;
            
            % Path speed (projection of velocity on path tangent)
            s_dot = x_current(4)*cos(psi_e) + x_current(5)*sin(psi_e);
            
            % Desired yaw rate
            r_des = delta_dot - beta_dot - obj.k_r*(psi_e - delt) + c_c*s_dot;
            
            % 7. Construct reference state
            x_ref = [lookahead_pt(1);       % Target x position
                    lookahead_pt(2);       % Target y position
                    psi_path + delt;      % Target heading (path + LOS)
                    0;          % Maintain surge velocity
                    0.1 * sign(y1);       % Small sway for correction 
                    r_des];               % Desired yaw rate
            
            % 9. Solve MPC
            [u_opt, ~] = obj.solve_mpc(x_current, x_ref);
            
            % 10. Store and visualize
            % obj.visualize(x_current, lookahead_pt, u_opt, y1);

            % Distribute thrust to vertical thrusters
            Tht = pi/4; %--> Thruster mounting angle
            j=[ cos(Tht) , cos(Tht), cos(Tht) , cos(Tht);
                -sin(Tht), sin(Tht), -sin(Tht), sin(Tht);
                -0.375   , 0.375   , 0.375    , -0.375];
         
            j_inv = pinv(j);

            % Apply thrust limits
            f_hor = j_inv * u_opt;
            f_hor = min(max(f_hor, -obj.Tlim), obj.Tlim);
        end
   
        function [lookahead_pt, s] = find_lookahead_point(obj, x_current)
            % Binary search for point exactly 'lookahead' distance ahead
            s_min = obj.current_s;
            s_max = min(obj.current_s + 0.5, obj.max_s); % Search 0.5m ahead max
            
            target_dist = obj.lookahead;
            tol = 0.01;
            
            while (s_max - s_min) > tol
                s_mid = (s_min + s_max)/2;
                pt = ppval(obj.path_spline, s_mid);
                dist = sqrt((pt(1)-x_current(1))^2+(pt(2)-x_current(2))^2);
                
                if dist < target_dist
                    s_min = s_mid;
                else
                    s_max = s_mid;
                end
            end
            
            s = (s_min + s_max)/2;
            lookahead_pt = ppval(obj.path_spline, s);
        end

        function [closest_pt, s] = find_closest_point(obj, x_current)
            % Find closest point on path to current position
            options = optimset('TolX', 1e-4);
            s = fminbnd(@(s) sqrt(sum((ppval(obj.path_spline, s) - x_current(1:2)').^2)),...
                       obj.current_s,...
                       min(obj.max_s, obj.current_s+0.2),...
                       options);
            closest_pt = ppval(obj.path_spline, s);
        end

        function visualize(obj, x, ref_pt, u, y1)
                figure(obj.figure_handle);
                hold on;
                
                % Clear previous dynamic elements
                h = findobj(gca, 'Tag', 'dynamic');
                if ~isempty(h), delete(h); end
                
                % Plot current state
                plot(x(1), x(2), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'Tag', 'dynamic');
                quiver(x(1), x(2), 0.5*cos(x(3)), 0.5*sin(x(3)), 'g', 'LineWidth', 2, 'Tag', 'dynamic');
                
                % Plot reference point
                plot(ref_pt(1), ref_pt(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'Tag', 'dynamic');
                
                % Display info
                text(0.05, 0.9, sprintf('Cross-track error: %.2f m', y1), ...
                    'Units', 'normalized', 'Tag', 'dynamic');
                text(0.05, 0.85, sprintf('Controls: Tx=%.1f N, Ty=%.1f N, Mz=%.1f Nm', ...
                    u(1), u(2), u(3)), 'Units', 'normalized', 'Tag', 'dynamic');
                
                drawnow;
            end

        function obj = update_visualization(obj, current_position)
            % Update the existing red marker (plot3 handle) with new AUV position
        
            if isvalid(obj.auv_marker)
                set(obj.auv_marker, ...
                    'XData', current_position(1), ...
                    'YData', current_position(2), ...
                    'ZData', current_position(3));
            else
                warning('AUV marker handle not found or invalid. Did you run set_path first?');
            end
        
            drawnow limitrate;  % Efficient redrawing
        end


        function [obj, pwm] = Actuate(obj, sens, n_des)
            %{ 
              Return Control pwms for eight thrusters
              """
              Inputs:
                    obj: class instance
                        class object
                    sens: structure
                        structure containing sensor readings (ORBSLAM + imu + pressure sensor)
                    n_des : double array
                        desired pose
              """
              Outputs:
                    obj: class instance
                        class object
                    pwm: 1x8 int array
                        pwm values sent to ESCs (1100 - 1900)
               """
            %}
            % Stabilization
            x=[sens.ps;sens.imu_euler(1);sens.imu_euler(2)];

            % Tracking
            des=[n_des(1), n_des(2), n_des(4),0,0,0];
            s=[sens.orbX, sens.orbY, sens.imu_euler(3)];
            [obj,s_dot]=obj.exact_differentiator(s,obj.s_prev);
            obj.s_prev = s;
            if isempty(obj.path_points)
                [obj, f_hor] = Track(obj, des, [s,s_dot]);
                y=[n_des(3);0;0];
                [obj, f_ver] = Stabilize(obj,x ,y);
            else
                obj = obj.update_visualization([sens.orbX, sens.orbY, sens.ps]);
                [obj, f_hor] = track_path(obj,[s,s_dot]);
                y=[obj.path_Z(obj.idx);0;0];
                [obj,f_ver] = Stabilize(obj,x ,y);
            end
            f_ver = reshape(f_ver,[1,4]);
            f_hor = reshape(f_hor,[1,4]);

            % Update states and controls
            obj.n = [sens.orbX, sens.orbY, sens.ps, sens.imu_euler(1), sens.imu_euler(2), ...
                sens.imu_euler(3)];
            obj.taw = [f_hor, f_ver];

            % PWM production
            pwm = T200_12v(obj.taw);
        end
    end
end

function pwm = T200_12v(F)
    %{ 
      Convert thrust to pwm according to T200 datasheet
      """
      Inputs:
            F: double
                thruster forces
      """
      Outputs:             
            pwm: int
                corresponding pwm values
       """
    %}
    pwm = zeros(size(F)); %--> Output pwm
    fit_poly_pos = [19.302480858079193, 39.3358541272889, -8.263830702456188, ...
        1.282961416635573, -0.112518356922797, 0.005724304685154, ...
        -1.678348654788274e-04, 2.629608127196461e-06, -1.703456825754830e-08]; 
            %--> Positive side fitting polynomial
    fit_poly_neg = [-19.25710704717055, 48.34221805064397, 12.322060404477401, ...
        2.354548958757111, 0.255975400767933, 0.016226987595385, ...
        5.950524882854001e-04, 1.169477093694967e-05, 9.526268034372109e-08]; 
            %--> Negative side fitting polynomial

    % Conversion formula
    for i = 1 : size(F, 2)
        if(F(i) > 0)
            pwm(i) = sum(fit_poly_pos .* F(i).^(0:8));
        elseif(F(i) < 0)
            pwm(i) = sum(fit_poly_neg .* F(i).^(0:8));
        else
            pwm(i) = 0;
        end
    end
    
    % Limit, shift and discretize pwm for ESCs
    pwm = round(min(max(pwm, -400), 400)) + 1500;
end