

clear
clc
close all

%% define parameters and construct a robot
link = [1;1];
com = [0.5;0];
mass = [1;1];
% inertia = 1/12*[mass(1)*link(1)^2;mass(2)*link(2)^2];
inertia = 1/12*[mass(1)*link(1)^2;0.5];
k = 0.3;
robot = UnderactuatedPlanarTwoLink(link,com,mass,inertia,k);

%% define simulation characteristics
t_end = 2;
ic = [pi/3 0 pi/4 0];
ut = transpose(0:0.05:t_end);
u = zeros(length(ut),1);
options = odeset('RelTol',1e-4,'AbsTol',1e-6*ones(1,4));
[T,X] = robot.Motion(ut,u,t_end,ic,options);

%% verify energy (for a free fall, the energy should conserve)
q = X(:,[1,3]);
q_dot = X(:,[2,4]);

E = zeros(length(T),1);
KE = zeros(length(T),1);
PE = zeros(length(T),1);
for i = 1:length(T)
    KE(i) = robot.kineticEnergy(q(i,:),q_dot(i,:));
    PE(i) = robot.potentialEnergy(q(i,:));
    E(i) = robot.totalEnergy(q(i,:),q_dot(i,:));
end

% plot energy
figure
plot(T,[E,KE,PE]);
legend('Total Energy','Kinetic Energy','Potential Energy')

%% animate motion
animationSwitch = 0;
if animationSwitch == 1
    figure
    for i = 1:length(T)
        robot = robot.setJointAngle(q(i,:));
        if i < length(T)
            dt = T(i+1)-T(i);
        end
        robot.animateMotion(dt)
        clf
    end
end

%% set up a motion planner
ts = 0; tf = 2;
t = ts:0.01:tf;
ic = [0 0 0 0];
fc = [pi/3 0 pi/4 0];
[yds,yds_dot,yds_ddot,yds_tdot,yds_qdot] = robot.motionPlanning(ic,fc,ts,tf);
[y,v] = robot.motionEvaluation(yds,yds_dot,yds_ddot,yds_tdot,yds_qdot,t);
u = robot.solveInputTorque(y,v);

options = odeset('RelTol',1e-4,'AbsTol',1e-6*ones(1,4));

[T1,X1] = robot.Motion(t,u,tf,ic,options);
q = X1(:,[1,3]);

animationSwitch = 1;
if animationSwitch == 1
    figure
    for i = 1:length(T1)
        robot.setJointAngle(q(i,:));
        if i < length(T1)
            dt = T1(i+1)-T1(i);
        end
        robot.animateMotion(dt)
        if i ~= length(T1)
            clf
        end
    end
end