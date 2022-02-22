classdef PlanarTwoLink < handle
    % PLANARTWOLINK defines a generic two-link planar manipulator. Each
    % link rotates about a revolute joint.
    % This program solves the kinematics and dynamics of system.
    %
    % Author:      Haohan Zhang
    % Edits:       Rand Hidayah
    % Affiliation: ROAR @ Columbia
    % Date:        5/1/2019
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ATTENTION SPRING 2020 Class
    %Don't do anythign to me.
    %All missing fucntions are in UnderActuatedPlanarTwoLink.m
    %It just inherits from this class.
    %Good luck!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    properties
        Link      % link lengths; a vector.
        COM       % location of center of mass; between 0-1; percentage of each full link length; a vector.
        Mass      % mass of each link; a vector.
        Inertia   % inertia of each link; a vector.
    end
    properties (Access = private)
        JointAngle    % angles of each joint, a vector
        JointTorque   % torque applied to each joint, a vector
    end
    
    methods
        % constructor
        function this = PlanarTwoLink(link,com,mass,inertia)
            % construct a robot
            if nargin > 0
                this.Link = link;
                this.COM = com;
                this.Mass = mass;
                this.Inertia = inertia;
            end
        end
        
        % setters
        function setJointAngle(this,value)
            % this assigns the joint angles
            this.JointAngle = value;
        end
        function setJointTorque(this,value)
            % this assigns the joint torques
            this.JointTorque = value;
        end
        
        % getters
        function value = getJointAngle(this)
            % this returns the joint angles
            value = this.JointAngle;
        end
        function value = getJointTorque(this)
            % this returns the joint torques
            value = this.JointTorque;
        end
        
        % helpers
        function A = calcPosA(this,q1)
            % this calculate the postion of the tip of the first link
            L = this.Link;
            x = L(1)*cos(q1);
            y = L(1)*sin(q1);
            A = [x;y];
        end
        function B = calcPosB(this,q)
            % this calculate the position of the far end of the 2nd link
            L = this.Link;
            x = L(1)*cos(q(1)) + L(2)*cos(q(1)+q(2));
            y = L(1)*sin(q(1)) + L(2)*sin(q(1)+q(2));
            B = [x;y];
        end
        
        % ode
        function [T,X] = Motion(this,t_end,ic,options)
            % this integrate the motion according to equations of motion
            [T,X] = ode45(@(t,x) equationOfMotion(t,x,this),[0 t_end],ic,options);
        end
        function dx = equationOfMotion(t,x,this)
            % this specifies the equations of motion of the system (ODE)
            % equation of motion:
            % A*q_ddot + B*q_dot + C = 0
            % convert into state space.
            % x1 = q1
            % x2 = q1_dot (= y1_dot)
            % x3 = q2
            % x4 = q2_dot (= y3_dot)
            % x_dot = D*y + E
            
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
            u = this.getJointTorque;
             
            % construct matries
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
             
            C1 = m1*g*L_C1*cos(x(1)) + m2*g*L1*cos(x(1)) + m2*g*L_C2*cos(x(1) + x(3)) - u(1);
            C2 = m2*g*L_C2*cos(x(1) + x(3)) - u(2);
            C = [C1;C2];
            
            D = -A_inv*B; 
            E = -A_inv*C;
            
            % ODE
            dx = zeros(4,1);
            dx(1) = x(2);
            dx(2) = D(1,1)*x(2) + D(1,2)*x(4) + E(1);
            dx(3) = x(4);
            dx(4) = D(2,1)*x(2) + D(2,2)*x(4) + E(2);
        end
        
        % energy
        function K = kineticEnergy(this,q,q_dot)
            % this computes the kinetic energy of the system
           
            % parameters
            M = this.Mass;
            I = this.Inertia;
            L = this.Link;
            C = this.COM;
            m1 = M(1); m2 = M(2);
            I_C1 = I(1); I_C2 = I(2);
            L1 = L(1); L2 = L(2);
            L_C1 = L1*C(1);L_C2 = L2*C(2);
            
            % calculation
            KE1 = 1/2*(m1*L_C1^2 + I_C1)*q_dot(1)^2;
            q12_dot = q_dot(1) + q_dot(2);
            KE2 = 1/2*m2*(L1^2*q_dot(1)^2 + L_C2^2*q12_dot^2 + ...
                2*L1*L_C2*cos(q(2))*q_dot(1)*q12_dot) + 1/2*I_C2*q12_dot^2;
            K = KE1 + KE2;
        end
        function V = potentialEnergy(this,q)
            % this returns the potential energy
            
            g = 9.8;
            % parameters
            M = this.Mass;
            I = this.Inertia;
            L = this.Link;
            C = this.COM;
            m1 = M(1); m2 = M(2);
            L1 = L(1); L2 = L(2);
            L_C1 = L1*C(1);L_C2 = L2*C(2);
            
            % calculation
            V = m1*g*L_C1*sin(q(1)) + m2*g*(L1*sin(q(1)) + L_C2*sin(q(1)+q(2)));
        end
        function E = totalEnergy(this,q,q_dot)
            % this returns total energy
            KE = this.kineticEnergy(q,q_dot);
            PE = this.potentialEnergy(q);
            E = KE + PE;
        end
        
        % visual
        function plotRobot(this)
            % this plot the geometry of robot
            O = [0;0];
            q = this.getJointAngle;
            A = this.calcPosA(q(1));
            B = this.calcPosB(q);
            L1 = [O A].'; L2 = [A B].';
            hold on
            line(L1(:,1),L1(:,2),'linewidth',2,'color','b');
            line(L2(:,1),L2(:,2),'linewidth',2,'color','b');
            plot(O(1),O(2),'k.','markersize',25);
            plot(A(1),A(2),'k.','markersize',25);
            hold off
        end
        function animateMotion(this,dt)
            % this animates the motion of the robot
            axis([-3 3 -3 3]);
            this.plotRobot;
            drawnow
            pause(dt) % pause with a 'correct' timing
%             clf
        end
    end
    
end

