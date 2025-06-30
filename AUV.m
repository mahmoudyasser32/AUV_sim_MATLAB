%{ 
    //Powered by ASmarine//
    AUV dynamic model class:
        This class can be used to create instances of an eight-thruster 
        configuration AUV with all hydrostatic and hydrodynamic parameters; 
        as well as sensors and actuators.
    Author: Samer A. Mohamed (https://github.com/SamMans) 
%}
classdef AUV
    properties(SetAccess = private)
    % Constants
        % Hydrostatic constants
            m %--> AUV mass (Unit: kg)
            V %--> AUV volume (Unit: m3)
            I0 %--> AUV inertia tensor (3x3 matrix, Unit: kg.m2)
            COB %--> Vector from COM to COB (1x3 vector, Unit: m)  
        % Hydrodynamic constants
            CD %--> Drag coefficients matrix (12x6 matrix, Unit: N/A)
            MA %--> Added mass coefficients matrix (12x6 matrix, Unit: N/A)
            CF %--> Lift coefficients matrix (12x6 matrix, Unit: N/A)
        % Thruster constants
            th_P2T %--> PWM to thrust mapping equation coefficients
            th_proj %--> Thruster projection matrix (6x8 matrix, Unit: m)
        % Depth sensor constants
            ps_res %--> Pressure sensor resolution (Unit: m)
            ps_relacc %--> Pressure sensor relative accuracy (Unit: N/A)
        
    % Variables
        % Kinematic variables
            n %--> Global pose vector (1x6 vector, Units: m, rad)
            v %--> Local velocity vector (1x6 vector, Units: m/s, rad/s)
        % Clock variables
            time %--> Simulation time (Unit: seconds)
            dT %--> Adaptive sampling time (Unit: seconds)
            dT_con %--> Controller sampling time (Unit: seconds)
    end
    methods
        function obj = AUV(model, ni, dt_c) 
            %{ 
              Constructor: Initializes AUV constants and initial states
              """
              Inputs:
                    model: string
                        AUV model name
                    ni: 1x6 array 
                        Initial global pose in meters - rads
                    dt_c : double
                        Controller sampling time
              """
              Outputs:
                    obj: class object 
                        Instance of AUV class 
               """
            %}
            % Read AUV parameters from .mat file
            if ~exist(pwd + "\Libraries\AUVs\" + model, 'dir')
               % Throw error if user chooses an unavailable model
               error(model + " doesn't exist, check the available models on SamSim libraries");
            end
            load(pwd + "\Libraries\AUVs\" + model + "\physical.mat");

            % Load Hydrostatic parameters
            obj.m = m_AUV; obj.V = V_AUV;
            obj.I0 = I0_AUV; obj.COB = COB_AUV;

            % Load Hydrodynamic parameters
            obj.MA = MA_AUV; 
            obj.CD = CD_AUV; 
            obj.CF = CF_AUV;
            
            % Load thruster constants
            Tht = pi/4; %--> Thruster mounting angle
            obj.th_P2T = p2t;
            obj.th_proj = [cos(Tht), cos(Tht), cos(Tht), cos(Tht), 0, 0, 0, 0;
                        -sin(Tht), sin(Tht), -sin(Tht), sin(Tht), 0, 0, 0, 0;
                        0, 0, 0, 0, 1, 1, 1, 1;
                        0, 0, 0, 0, Ry_AUV, -Ry_AUV, -Ry_AUV, Ry_AUV;
                        0, 0, 0, 0, -Rx_AUV, -Rx_AUV, Rx_AUV, Rx_AUV;
                        -R_Hor_AUV, R_Hor_AUV, R_Hor_AUV, -R_Hor_AUV, 0, 0, 0, 0];
            
            % Load pressure sensor constants
            obj.ps_res = psres;
            obj.ps_relacc = psrelacc; 
                          
            % Initial conditions
            obj.n = ni; %--> global states
            obj.v = zeros(1, 6); %--> local states
            
            % Simulation time initialization
            obj.time = 0; obj.dT_con = dt_c;
            obj.dT = obj.dT_con;
        end
        function state = getstate(obj)
            %{ 
              Return current AUV states with associated time stamp
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
              """
              Outputs:
                    state: struct 
                        AUV states structure with three properties 
               """
            %}
            state = [obj.n, obj.v, obj.time];
        end
        function t = gettime(obj)
            %{ 
              Return current simulation time
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
              """
              Outputs:
                    t: double
                        Simulation time
               """
            %}
            t = obj.time;
        end
        function sens = getsens(obj)
            %{ 
              Return sensor readings
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
              """
              Outputs:
                    sens: structure
                        Sensor reading structure
               """
            %}
            sens.orbX = obj.n(1) ;%+ normrnd(0, 1/30); %--> Visual SLAM X reading
            sens.orbY = obj.n(2) ;%+ normrnd(0, 1/30); %--> Visual SLAM Y reading
            sens.ps = obj.pressure_sensor(); %--> Pressure sensor reading
            sens.imu_euler = obj.n(4:6); %--> IMU reading
        end
        function J = Jacobian(obj)
            %{ 
              Returns transformation jacobian from local to global frame
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
              """
              Outputs:
                    J: 6x6 matrix
                        AUV jacobian
               """
            %}
            J = [cos(obj.n(6))*cos(obj.n(5)) -sin(obj.n(6))*cos(obj.n(4))+cos(obj.n(6))*sin(obj.n(5))*sin(obj.n(4))...
                sin(obj.n(6))*sin(obj.n(4))+cos(obj.n(6))*cos(obj.n(4))*sin(obj.n(5)) 0 0 0;
                sin(obj.n(6))*cos(obj.n(5)) cos(obj.n(6))*cos(obj.n(5))+sin(obj.n(4))*sin(obj.n(5))*sin(obj.n(6)) ...
                -cos(obj.n(6))*sin(obj.n(4))+sin(obj.n(5))*sin(obj.n(6))*cos(obj.n(4)) 0 0 0;
                -sin(obj.n(5)) cos(obj.n(5))*sin(obj.n(4)) cos(obj.n(5))*cos(obj.n(4)) 0 0 0;
                0 0 0 1 sin(obj.n(4))*tan(obj.n(5)) cos(obj.n(4))*tan(obj.n(5))
                0 0 0 0 cos(obj.n(4)) -sin(obj.n(4));
                0 0 0 0 sin(obj.n(4))/cos(obj.n(5)) cos(obj.n(4))/cos(obj.n(5))];
        end
        function Th = convert_pwm(obj, PWM)
            %{ 
              Convert PWM signals to thrust forces
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
                    PWM: 1x6 array 
                        PWM signals
              """
              Outputs:
                    Th: 1x6 array
                        Thrust forces
              """
            %}
            PWM = round(PWM) - 1500;
            PWM = min(max(PWM, -400), 400); %--> PWM range limiter
            Th = zeros(1, 8);
            for i = 1 : size(obj.th_P2T, 2) / 2
                Th = Th + ((PWM >= 0) * obj.th_P2T(i) + ...
                    (1 - (PWM >= 0)) * obj.th_P2T(size(obj.th_P2T, 2)/2 + i)) .* PWM.^(i - 1);
            end
            Th = Th .* (PWM ~= 0); %--> Make sure zero pwms are mapped to zero thrust
        end
        function [vdot, ndot] = Newton(obj, PWM)
            %{ 
              Compute vehicle accelerations using Newton's 2nd law
              """
              Inputs:
                    PWM: 1x6 array 
                        PWM signals
                    obj: class object 
                        Instance of AUV class
              """
              Outputs:
                    vdot: 1x6 array
                        Local acceleration
                    ndot: 1x6 array
                        Global pose rate of change
              """
            %}
            % Calculate Hydrostatic forces
            g = 9.81; raw = 1000;
            Tg=  [-sin(obj.n(5)) * obj.m * g;
                  cos(obj.n(5)) * sin(obj.n(4)) * obj.m * g;
                  cos(obj.n(5)) * cos(obj.n(4)) * obj.m * g;
                  0;
                  0;
                  0]; %--> Weight force/moment vector w.r.t local frame         
            Tb=  [sin(obj.n(5)) * raw * obj.V * g;
                  -cos(obj.n(5)) * sin(obj.n(4)) * raw * obj.V * g;
                  -cos(obj.n(5)) * cos(obj.n(4)) * raw * obj.V * g;
                  -obj.COB(2) * cos(obj.n(5)) * cos(obj.n(4)) * raw * obj.V * g + obj.COB(3) * cos(obj.n(5)) * sin(obj.n(4)) * raw * obj.V * g;
                  obj.COB(3) * sin(obj.n(5)) * raw * obj.V * g + obj.COB(1) * cos(obj.n(5)) * cos(obj.n(4)) * raw * obj.V * g;
                  -obj.COB(1) * cos(obj.n(5)) * sin(obj.n(4)) * raw * obj.V * g - obj.COB(2) * sin(obj.n(5)) * raw * obj.V * g]; 
                      %--> Buoyancy force/moment vector w.r.t local frame 
            T_Static = Tg + Tb; %--> Total Hydrostatic Force/Moment w.r.t local frame 

            % Determine current drag, lift and added mass matrices 
            dir = zeros(6,6); 
            for i = 1:6
                dir(i,:) = obj.v >= 0; %--> Motion direction (+ve or -ve)
            end
            cd = dir.*obj.CD(1:6,:) + (1 - dir).*obj.CD(7:12,:);
            cf = dir.*obj.CF(1:6,:) + (1 - dir).*obj.CF(7:12,:);
            ma = dir.*obj.MA(1:6,:) + (1 - dir).*obj.MA(7:12,:);
            
            % Calculate Hydrodynamic forces
            Td = -cd * (obj.v.'.*abs(obj.v.')); %--> Drag force/moment vector w.r.t local frame    
            Tl = cf * (obj.v.'.*obj.v.'); %--> Lift force/moment vector w.r.t local frame 
            T_Dynamic = Td + Tl; %--> Total Hydrodynamic Force/Moment w.r.t local frame 
            
            % Calculate Thrust forces
            Tt = obj.th_proj * convert_pwm(obj, PWM).'; %--> Resultant thrust force/moment w.r.t local frame 
            
            % Newton's 2nd law
            M = [diag([obj.m, obj.m, obj.m]), zeros(3, 3);
               zeros(3, 3), obj.I0];  %--> System inertia tensor

            Crb = [zeros(3,3), (-obj.m * skew(obj.v(1:3).'));  
                   (-obj.m * skew(obj.v(1:3).')), (-skew(obj.I0 * obj.v(4:6).'))]; %--> Coriolis

            Ca = [zeros(3,3), (-skew(ma(1:3,1:3) * obj.v(1:3).' + ma(1:3,4:6) * obj.v(4:6).'));  
                 (-skew(ma(1:3,1:3) * obj.v(1:3).' + ma(1:3,4:6) * obj.v(4:6).')),...
                 (-skew(ma(4:6,1:3) * obj.v(1:3).' + ma(4:6,4:6) * obj.v(4:6).'))]; %--> Added Mass Coriolos
            T = T_Static + T_Dynamic + Tt; %--> Total acting forces/moments
            D=0;
            % if obj.time >= 20 && obj.time <= 100
                D=10*sin(2*pi*2*obj.time); %disturbance
            % end
            vdot = ((M + ma) \ (T + D - (Ca + Crb) * obj.v.')).'; %--> AUV acceleration
            ndot = (Jacobian(obj) * obj.v.').'; %--> Derivative of global states
        end
        function [obj] = Solve(obj, PWM)
            %{ 
              Solves the dynamic model and computes state evolution
              using fourth-order Runge-Kutta method
              """
              Inputs:
                    obj: class object 
                        Instance of AUV class
                    PWM: 1x6 array 
                        PWM signals
              """
              Outputs:
                    obj: class object 
                        Instance of AUV class
              """
            %}
            % Ranga-kutta parameters 
            dn_max = 0.01; %--> max relative change tolerance
            dn_min = 0.001; %--> min relative change tolerance
            dT_min = 10^-4; %--> Min possible simulation time step
            
            % Save starting states 
            n_start = obj.n; %--> Starting global pose
            v_start = obj.v; %--> Starting local velocity
            
            % Normal step (apply Runge-Kutta method for normal time step)
            [vdot_1, ndot_1] = Newton(obj, PWM); %--> Get acceleration (K1 slope)
            avgdn = ndot_1; avgdv = vdot_1; %--> Averaging of slopes (currently add K1 slope only)
            obj.n = n_start + (obj.dT / 2) * ndot_1; %--> Partial pose evolution 
            obj.v = v_start + (obj.dT / 2) * vdot_1; %--> Partial velocity evolution
            [vdot_2, ndot_2] = Newton(obj, PWM); %--> Get acceleration (K2 slope)
            avgdn = avgdn + 2 * ndot_2; avgdv = avgdv + 2 * vdot_2; %--> Averaging of slopes (currently K1 & K2 slopes)
            obj.n = n_start + (obj.dT / 2) * ndot_2; %--> Partial pose evolution 
            obj.v = v_start + (obj.dT / 2) * vdot_2; %--> Partial velocity evolution
            [vdot_3, ndot_3] = Newton(obj, PWM); %--> Get acceleration (K3 slope)
            avgdn = avgdn + 2 * ndot_3; avgdv = avgdv + 2 * vdot_3; %--> Averaging of slopes (currently for K1, K2 & K3 slopes)
            obj.n = n_start + obj.dT * ndot_3; %--> Partial pose evolution 
            obj.v = v_start + obj.dT * vdot_3; %--> Partial velocity evolution
            [vdot_4, ndot_4] = Newton(obj, PWM); %--> Get acceleration (K4 slope)
            avgdn = avgdn + ndot_4; avgdv = avgdv + vdot_4; %--> Averaging of slopes (currently for all slopes)
            step_n = n_start + (obj.dT / 6) * avgdn; %--> Compute normal step global pose update
            step_v = v_start + (obj.dT / 6) * avgdv; %--> Compute normal step local velocity update
            obj.n = n_start; obj.v = v_start; %--> Reset states to their initial values
            
            % Half step (apply Runge-Kutta method for half time step)
            [vdot_1, ndot_1] = Newton(obj, PWM); %--> Get acceleration (K1 slope)
            avgdn = ndot_1; avgdv = vdot_1; %--> Averaging of slopes (currently add K1 slope only)
            obj.n = n_start + (obj.dT / 4) * ndot_1; %--> Partial pose evolution 
            obj.v = v_start + (obj.dT / 4) * vdot_1; %--> Partial velocity evolution
            [vdot_2, ndot_2] = Newton(obj, PWM); %--> Get acceleration (K2 slope)
            avgdn = avgdn + 2 * ndot_2; avgdv = avgdv + 2 * vdot_2; %--> Averaging of slopes (currently K1 & K2 slopes)
            obj.n = n_start + (obj.dT / 4) * ndot_2; %--> Partial pose evolution 
            obj.v = v_start + (obj.dT / 4) * vdot_2; %--> Partial velocity evolution
            [vdot_3, ndot_3] = Newton(obj, PWM); %--> Get acceleration (K3 slope)
            avgdn = avgdn + 2 * ndot_3; avgdv = avgdv + 2 * vdot_3; %--> Averaging of slopes (currently for K1, K2 & K3 slopes)
            obj.n = n_start + (obj.dT / 2) * ndot_3; %--> Partial pose evolution
            obj.v = v_start + (obj.dT / 2) * vdot_3; %--> Partial velocity evolution
            [vdot_4, ndot_4] = Newton(obj, PWM); %--> Get acceleration (K4 slope)
            avgdn = avgdn + ndot_4; avgdv = avgdv + vdot_4; %--> Averaging of slopes (currently for all slopes)
            half_step_n = n_start + (obj.dT / 12) * avgdn; %--> Compute half step global pose update
            half_step_v = v_start + (obj.dT / 12) * avgdv; %--> Compute half step local velocity update
            obj.n = n_start; obj.v = v_start; %--> Reset states to their initial values

            % Double step (apply Runge-Kutta method for double time step)
            [vdot_1, ndot_1] = Newton(obj, PWM); %--> Get acceleration (K1 slope)
            avgdn = ndot_1; avgdv = vdot_1; %--> Averaging of slopes (currently add K1 slope only)
            obj.n = n_start + obj.dT * ndot_1; %--> Partial pose evolution 
            obj.v = v_start + obj.dT * vdot_1; %--> Partial velocity evolution
            [vdot_2, ndot_2] = Newton(obj, PWM); %--> Get acceleration (K2 slope)
            avgdn = avgdn + 2 * ndot_2; avgdv = avgdv + 2 * vdot_2; %--> Averaging of slopes (currently K1 & K2 slopes)
            obj.n = n_start + obj.dT * ndot_2; %--> Partial pose evolution
            obj.v = v_start + obj.dT * vdot_2; %--> Partial velocity evolution
            [vdot_3, ndot_3] = Newton(obj, PWM); %--> Get acceleration (K3 slope)
            avgdn = avgdn + 2 * ndot_3; avgdv = avgdv + 2 * vdot_3; %--> Averaging of slopes (currently for K1, K2 & K3 slopes)
            obj.n = n_start + (obj.dT * 2) * ndot_3; %--> Partial pose evolution
            obj.v = v_start + (obj.dT * 2) * vdot_3; %--> Partial velocity evolution
            [vdot_4, ndot_4] = Newton(obj, PWM); %--> Get acceleration (K4 slope)
            avgdn = avgdn + ndot_4; avgdv = avgdv + vdot_4;  %--> Averaging of slopes (currently for all slopes)
            dbl_step_n = n_start + (obj.dT / 3) * avgdn; %--> Compute double step global pose update
            dbl_step_v = v_start + (obj.dT / 3) * avgdv; %--> Compute double step local velocity update
            
            % Choose correct step size, update states & update simulation time
            if(obj.dT < dT_min) %--> Check sampling time doesn't fall below a minimum threshold
                % If it does, stick to normal step state evolution
                % reset sampling time to minimum sampling time
                obj.n = step_n; obj.v = step_v;
                obj.time = obj.time + obj.dT; %--> Update simulation time
                obj.dT = dT_min;
            else
                if(any(abs(step_n - half_step_n) > dn_max)) %--> Use half step if step/half step error exceeds threshold "dn_max"
                    obj.dT = obj.dT / 2;
                    obj.n = half_step_n; obj.v = half_step_v;
                    obj.time = obj.time + obj.dT; %--> Update simulation time
                else
                    if(all(abs(step_n - dbl_step_n) < dn_min) && obj.dT_con >= obj.dT * 2) %--> Use dbl step if error step/double step is less than threshold (dn_min)
                        obj.dT = obj.dT * 2;
                        obj.n = dbl_step_n; obj.v = dbl_step_v;
                        obj.time = obj.time + obj.dT; %--> Update simulation time
                    else
                        obj.n = step_n; obj.v = step_v; %--> Stick to normal step if none
                        obj.time = obj.time + obj.dT; %--> Update simulation time
                    end
                end
            end
        end
        function noise_depth = pressure_sensor(obj)
            % PRESSURE_SENSOR Simulates the depth reading from a pressure sensor with noise.
            % 
            % INPUTS:
            %   depth - The true depth of the sensor (cm).
            %   pitch - The pitch angle of the sensor (radians).
            %
            % OUTPUT:
            %   noise_depth - The measured depth with added sensor noise (cm).
            
            dist = 25; % Distance between the center of gravity and the pressure sensor (cm)
            
            % Compute the ideal (noiseless) pressure sensor reading
            ideal_reading = obj.n(3)*100 + dist * sin(obj.n(5));
            
            % Generate random noise (uniform distribution between 0 and 0.12)
            noise = 0.12 * (-1 + 2 * rand); 
            
            % Add noise to the ideal reading to simulate sensor inaccuracy
            noise_depth = ideal_reading + noise;
            noise_depth = noise_depth/100;
        end
    end
end
% Ordinary functions
function S = skew(a)
    % Return a skew symmetric matrix made up of elements of vector a
    S = [0 -a(3,1) a(2,1); 
        a(3,1) 0 -a(1,1); 
        -a(2,1) a(1,1) 0];
end