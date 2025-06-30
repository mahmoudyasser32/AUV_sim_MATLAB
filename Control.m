classdef Control
    properties(SetAccess = private)
        % AUV parameters
        dT %--> Controller sampling time (1*1 double, unit: sec.)
        Tlim %--> Thruster output limit (1*1 double, unit: N)
        rh %--> AUV horizontal thruster arm length (1*1 double, unit: m)
        I_lin %--> AUV linear inertia
        D_lin %--> AUV linear drag
        I_rot %--> AUV rotational inertia
        D_rot %--> AUV rotational drag

        % Dynamic properties
        taw %--> Control thrust (1*8 array of double, unit: N)
        n %--> Global pose feedback (1*6 array of double, unit: m and/or rad)
        v %--> Local velocity estimate (1*6 array of double, unit: m/s and/or rad/s)

        % Stabilizer parameters
        Ksp %--> Proportional gain (1*1 double, unitless)
        Ksd %--> Derivative gain (1*1 double, unitless)

        % Trajectory tracking parameters
        solver %--> Non-linear MPC solver object
        N %--> MPC prediction horizon (1*1 double, unitless)
    end
    methods
        function obj = Control(sens, dT)
            %{ 
              Constructor
              """
              Inputs:
                    sens: structure
                        structure containing sensor readings (EKF/ORBSLAM + imu + pressure sensor)
                    dT : double
                        controller sampling time
              """
              Outputs:
                    obj: class instance
                        class object
               """
            %}
            % Initialize states and controls
            obj.taw = zeros(1,8);
            obj.n = [sens.orbX, sens.orbY, sens.ps, sens.imu_euler(1), sens.imu_euler(2), ...
                sens.imu_euler(3)]; %--> Initialize pose using raw sensor data
            obj.v = zeros(1,6); %--> Initialize local velocity estimate
            obj.dT = dT; obj.Tlim = 6; obj.rh = 0.375;

            % Initialize stabilizer parameters
            obj.Ksp = 1000; %--> Proportional gain
            obj.Ksd = 200; %--> Derivative gain

            % Initialize trajectory tracker parameters
            addpath(pwd + "\Dependencies\Casadi")
            import casadi.*
            obj.N = 20; %--> Instantaneous prediction horizon
            obj.I_lin = 25.305 + 360.29; %--> Linear inertia (mass + added mass)
            obj.D_lin = 298.28; %--> Linear drag
            obj.I_rot = 8.544; %--> Rotational inertia (moment of inertia + added mass)
            obj.D_rot = 2.8145; %--> Rotational drag

            % Nonlinear programming problem definition
            x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi'); %--> Global pose states
            u = SX.sym('u'); r = SX.sym('r'); %--> Local velocity states
            states = [x; y; psi; u; r]; %--> Model internal states vector
            n_states = length(states); %--> Number of internal states
            F = SX.sym('F'); T = SX.sym('T'); %--> Kinetic controls (Force and Torque)
            controls = [F; T]; %--> Model controls
            n_controls = length(controls); %--> Number of controls
            rhs = [u*cos(psi); u*sin(psi); r; (F - obj.D_lin*u*abs(u))/obj.I_lin; ...
                (T - obj.D_rot*r*abs(r))/obj.I_rot]; %--> States rate of change
            f = Function('f', {states,controls}, {rhs}); %--> Model function
            U = SX.sym('U', n_controls, obj.N); %--> Sequence of predictive controls
            P = SX.sym('P', n_states + 3); %--> Initial + desired states vector
            X = SX.sym('X', n_states, obj.N+1); %--> Sequence of state predictions
            X(:,1) = P(1:n_states); %--> Initial kinematic states
            for k = 1:obj.N
                st = X(:,k); con = U(:,k);
                f_value = f(st,con);
                st_next = st + (obj.dT * f_value); %--> State evolution
                X(:,k+1) = st_next;
            end
            loss = 0; %--> Objective function initialized to zero
            g = []; %--> Constraints on states (none)
            Q = zeros(3,3); Q(1,1) = 700000; Q(2,2) = 500000; Q(3,3) = 50000; %--> Performance weight matrix
            for k = 1 : obj.N
                st = X(1:3,k);
                loss = loss + (st-P(n_states+1:end))'*Q*(st-P(n_states+1:end)); %--> Add individual losses to objective function
            end
            OPT_variables = reshape(U, n_controls*obj.N, 1); %--> Optimization variables
            nlp = struct('x', OPT_variables, 'f', loss, 'g', g, 'p', P); %--> Nonlinear programming problem
            opts = struct; %--> Options structure
            opts.ipopt.max_iter = 100;
            opts.ipopt.print_level = 0;
            opts.print_time = 0;
            opts.ipopt.acceptable_tol = 1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
            obj.solver = nlpsol('solver','ipopt',nlp,opts); %--> Solver object
        end
        function f_ver = Stabilize(obj, Z_des, Z_curr)
            %{ 
              Return stabilization thrusts for 4 vertical thrusters
              """
              Inputs:
                    obj: class instance
                        class object
                    Z_des: double
                        desired depth
                    Z_curr : double
                        current depth reading
              """
              Outputs:             
                    f_ver: 1x4 double array
                        vertical thrust values
               """
            %}
            % Depth error calculation
            e = Z_des - Z_curr; %--> Error

            % Discretized control law
            e_old = Z_des - obj.n(3);
            F_old = sum(obj.taw(5:8));
            F = (obj.Ksp + (2*obj.Ksd)/obj.dT) * e + (obj.Ksp - (2*obj.Ksd)/obj.dT) * e_old - F_old; %--> Control force
            f = F/4; %--> Force per thruster (four thrusters, so divide by 4)
            f = min(max(f, -obj.Tlim), obj.Tlim); %--> Thrust limiting
            f_ver = f * ones(1, 4);
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
            % Set nonlinear programming problem arguments
            args = struct; Fmax = 2*cos(pi/4)*obj.Tlim; Tmax = 2*obj.rh*obj.Tlim;
            args.lbx(1:2:2*obj.N-1,1) = -Fmax; args.lbx(2:2:2*obj.N,1) = -Tmax;
            args.ubx(1:2:2*obj.N-1,1) = Fmax; args.ubx(2:2:2*obj.N,1) = Tmax;
            U0 = zeros(obj.N, 2);  %--> Predictive control inputs 
            args.p = [pose_curr.'; obj.v(1); obj.v(6); pose_des.']; %--> Set the values of the parameters vector
            args.initialStates = reshape(U0', 2*obj.N, 1); %--> Initial value of the optimization variables

            % Solve the nonlinear programming problem
            sol = obj.solver('x0', args.initialStates, 'lbx', args.lbx, 'ubx', args.ubx,...
                'p',args.p);
            U = reshape(full(sol.x)', 2, obj.N)';
            F = U(1,1); T = U(1,2); %--> Current kinetic controls
            obj.v(1) = obj.v(1) + ((F-obj.D_lin*obj.v(1)*abs(obj.v(1)))/obj.I_lin)*obj.dT;
            obj.v(6) = obj.v(6) + ((T-obj.D_rot*obj.v(6)*abs(obj.v(6)))/obj.I_rot)*obj.dT;

            % Calculate kinetic controls
            T1 = F/(2*cos(pi/4)); T2 = F/(2*cos(pi/4)); %--> Translational motion part
            T3 = T/(2*obj.rh); T4 = -T/(2*obj.rh); %--> Rotational motion part
            T1 = min(max(T1, -obj.Tlim), obj.Tlim); T2 = min(max(T2, -obj.Tlim), obj.Tlim); 
            T3 = min(max(T3, -obj.Tlim), obj.Tlim); T4 = min(max(T4, -obj.Tlim), obj.Tlim); %--> Thrust limiting
            f_hor = [T1, T2, T3, T4];
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
            f_ver = Stabilize(obj, n_des(3), sens.ps);

            % Tracking
            % [obj, f_hor] = Track(obj, [n_des(1), n_des(2), n_des(4)], ...
            %     [sens.orbX, sens.orbY, sens.imu_euler(3)]);
            f_hor = [0,0,0,0];

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