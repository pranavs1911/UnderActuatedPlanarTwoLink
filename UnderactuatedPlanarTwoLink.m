classdef UnderactuatedPlanarTwoLink < PlanarTwoLink
    %UNDERACTUATEDPLANARTWOLINK is a subclass of the generic planar two
    %link robot. Add more functionalities in this class so that it can
    % be used to plan motions for a robot with only a single actuator which
    % is mounted on the proximal joint. Define its parameters such that the
    % center of mass of the 2nd link is located at the end of the 1st link.
    % Add a torsional spring (spring constant k, zero free length) at the
    % second joint so that the coupling effect of the motion can be
    % chracterized.
    %
    % Author:       Rand Hidayah
    % Affiliation: ROAR @ Columbia
    % Date:        4/14/2020
    
    properties
        k   % spring constant of the torsional spring
    end
    
    methods
        
        % constructor
        function this = UnderactuatedPlanarTwoLink(link,com,mass,inertia,springConst)
            if nargin == 0
                link = [1;1];
                com = [0.5;0];
                mass = [1;1];
                inertia = 1/12*[mass(1)*link(1)^2;mass(2)*link(2)^2];
                springConst = 0;
            end
            this@PlanarTwoLink(link,com,mass,inertia);
            this.k = springConst;
        end
        
        % new equation of motion
        function [T,X] = Motion(this,ut,u1,t_end,ic,options)
            % this integrate the motion according to equations of motion
            [T,X] = ode45(@(t,x) equationOfMotion(t,x,ut,u1,this),[0 t_end],ic,options);
        end
        function dx = equationOfMotion(t,x,ut,u1,this)
            % this specifies the equations of motion of the system (ODE)
            % Note: You can copy the original ODE in the generic class and
            % make some minor changes.
            g = 9.8;
            % parameters
            M = this.Mass;
            I = this.Inertia;
            L = this.Link;
            C = this.COM;
            m1 = M(1); m2 = M(2);
            I_C1 = I(1); I_C2 = I(2);
            L1 = L(1); L2 = L(2);
            L_C1 = L1*C(1);L_C2 = L2*C(2);
            %You can mess with the control law%
            u1 = interp1(ut,u1,t);
            
            %end%
            K = this.k;
            
            % construct matrices
            A11 = m1*L_C1^2 + I_C1 + m2*L1^2 + m2*L_C2^2 + I_C2 + 2*m2*L1*L_C2*cos(x(3));
            A12 = m2*L1*L_C2*cos(x(3)) + I_C2 + m2*L_C2^2;
            A21 = m2*L_C2^2 + I_C2 + m2*L1*L_C2*cos(x(3));
            A22 = m2*L_C2^2 + I_C2;
            A = [A11 A12;A21 A22];
            A_inv = inv(A);
            
            B11 = -2*m2*L1*L_C2*sin(x(3))*x(4);
            B12 = -m2*L1*L_C2*sin(x(3))*x(4);
            B21 = m2*L1*L_C2*sin(x(3))*x(2);
            B22 = 0;
            B = [B11 B12;B21 B22];
            
            %%%%%%%%%%%%%%%%%%%%%
            %ATTENTION
            %ATTENTION
            %PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
            % Hint: don't forget the spring, or other torques.
            C1 = m1*g*L_C1*cos(x(1)) +m2*g*L1*cos(x(1)) + m2*g*L_C2*cos(x(1)+x(3)) ; %TODO Fill in your potential/input terms.
            C2 = m2*g*L_C2*cos(x(1) + x(3)) + K*(x(3)); %^
            C = [C1;C2];
            
            % calculation.
            D = -A_inv*B; % Don't change
            E = -A_inv*C; % Don't change
            %Aqddot + Bqdot + C = tau
            %qddot = Ainv*tau + Dqdot + E
            del = A11;
            gamma = C1/cos(x(1));
            q1ddot = (u1 + K*x(3) - gamma)/(del - I_C2);
            q2ddot = (I_C2*u1 - I_C2*gamma*cos(x(1)) + K*del*x(3));
            % ODE
            dx = zeros(4,1);
            dx(1) = x(2); % TODO Fill in ODE expressions
            dx(2) = q1ddot; % 
            dx(3) = x(4); % 
            dx(4) = q2ddot; %                                                                                                                   
            %%%%%%%%%%%%%%%%%%%%%
            %END OF PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
        end
        
        % potential energy with springs
        function V = potentialEnergy(this,q)
            % this returns the potential energy. 
            % Note there is a torsion spring added.
            g = 9.8;
            % parameters
            M = this.Mass;
            %I = this.Inertia;
            L = this.Link;
            C = this.COM;
            m1 = M(1); m2 = M(2);
            L1 = L(1); %L2 = L(2);
            L_C1 = L1*C(1);
            %L_C2 = L2*C(2);
            K = this.k;
            
            %%%%%%%%%%%%%%%%%%%%%
            %ATTENTION
            %ATTENTION
            %PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
            % calculation
            Vg = m1*g*L_C1*sin(q(:,1)) + m2*g*L1*sin(q(:,1)); % Fill in your expressions for Potential Energy
            Vs = 0.5*K*(q(:,3)).^2; % ^
            V = Vs+Vg;
            %%%%%%%%%%%%%%%%%%%%%
            %END OF PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
        end

        % plan desired motion using feedback linearization
        function y = transformation(this,x)
            % this transforms the original states into a new state space,
            % where, x = [x1 x2 x3 x4] = [q1 q1_dot q2 q2_dot], 
            %        y = [y1 y2 y3 y4] = [y y_dot y_ddot y_tdot],
            I = this.Inertia;
            K = this.k;
            %%%%%%%%%%%%%%%%%%%%%
            %ATTENTION
            %ATTENTION
            %PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
            I_C2 = I(2);
            y1 =I_C2*(x(1) + x(3));% Fill in state space transofrmation;
            y2 =I_C2*(x(2) + x(4));%^
            y3 = -K*x(3);%^
            y4 = -K*x(4);%^
            %%%%%%%%%%%%%%%%%%%%%
            %END OF PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
            y = [y1 y2 y3 y4];
        end
        
        % trajectory generator
        function [yds,yds_dot,yds_ddot,yds_tdot,yds_qdot] = motionPlanning(this,ic_x,fc_x,~,tf)
            % this plans a motion using a 7th order polynomial based on
            % inital and final conditions with respect to the transformed
            % states y, however, the function takes the initial and final
            % conditions of the original states. It returns the polynomial
            % expression symbolically.
            ic_y = this.transformation(ic_x);
            fc_y = this.transformation(fc_x);
            I = this.Inertia;
            K = this.k;
            %%%%%%%%%%%%%%%%%%%%%
            %ATTENTION
            %ATTENTION
            %PART TO BE FILLED IN
            A = [0 0 0 0 0 0 0 1;tf^7 tf^6 tf^5 tf^4 tf^3 tf^2 tf^1 1;0 0 0 0 0 0 1 0;7*tf^6 6*tf^5 5*tf^4 4*tf^3 3*tf^2 2*tf 1 0; 0 0 0 0 0 2 0 0; 42*tf^5 30*tf^4 20*tf^3 12*tf^2 6*tf^1 2 0 0; 0 0 0 0 6 0 0 0;210*tf^4 120*tf^3 60*tf^2 24*tf 6 0 0 0]; % TODO Fill in your matrix.
            b = [ic_y(1);fc_y(1);ic_y(2);fc_y(2);ic_y(3);fc_y(3);ic_y(4);fc_y(4)]; % TODO Initial conditions Vector
            par = A\b; % Don't change 
            % symbolic expression
            syms T
            yds      = par(1)*T^7 + par(2)*T^6 + par(3)*T^5 + par(4)*T^4 + par(5)*T^3 + par(6)*T^2 + par(7)*T^1 + par(8);% TODO Fill in your symbolic expression,
            yds_dot  = 7*par(1)*T^6 + 6*par(2)*T^5 + 5*par(3)*T^4 + 4*par(4)*T^3 + 3*par(5)*T^2 + 2*par(6)*T^1 + par(7);% ^ based on a seventh order polynomial.
            yds_ddot = 42*par(1)*T^5 + 30*par(2)*T^4 + 20*par(3)*T^3 + 12*par(4)*T^2 + 6*par(5)*T^1 + 2*par(6);% ^
            yds_tdot = 210*par(1)*T^4 + 120*par(2)*T^3 + 60*par(3)*T^2 + 24*par(4)*T^1 + 6*par(5);% ^
            yds_qdot = 840*par(1)*T^3 + 360*par(2)*T^2 + 120*par(3)*T^1 + 24*par(4);% ^
            %%%%%%%%%%%%%%%%%%%%%
            %END OF PART TO BE FILLED IN
            %%%%%%%%%%%%%%%%%%%%%
        end
        function [y,v] = motionEvaluation(~,yds,yds_dot,yds_ddot,yds_tdot,yds_qdot,t)
            % this evaluates the symbolic expression numerically
            % based on the time vector t.
            syms T
            yd = zeros(length(t),1);
            yd_dot = zeros(length(t),1);
            yd_ddot = zeros(length(t),1);
            yd_tdot = zeros(length(t),1);
            yd_qdot = zeros(length(t),1);
            
            for i = 1:length(t)
                yd(i) = double(subs(yds,T,t(i)));
                yd_dot(i) = double(subs(yds_dot,T,t(i)));
                yd_ddot(i) = double(subs(yds_ddot,T,t(i)));
                yd_tdot(i) = double(subs(yds_tdot,T,t(i)));
                yd_qdot(i) = double(subs(yds_qdot,T,t(i)));
            end
            
            y = [yd yd_dot yd_ddot yd_tdot];
            v = yd_qdot;
        end
        
        % solve for input torque u
        function u = solveInputTorque(this,y,v)
            % this computes u from v.
            I = this.Inertia; m = this.Mass; l1 = this.Link(1); lc1 = l1*this.COM(1);
            K = this.k;
            a1 = I(1)+I(2)+m(1)*lc1^2+m(2)*l1^2; a2 = I(2);
            a3 = I(2); a4 = I(2);
            
            u = zeros(length(v),1);
            for i = 1:length(v)
                x1 = 1/I(2)*y(i,1) + 1/K*y(i,3);
                b1 = (m(1)*lc1+m(2)*l1)*9.8*cos(x1);
                u(i) = -(a2*a3-a1*a4)/(K*a3)*v(i) + b1 + a1/a3*y(i,3);
            end
        end
        
    end
end