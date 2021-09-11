%% ROBOTICS - Aristotle University of Thessaloniki 
% Project - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Project Description:
% In this project the kinematics and trajectory planning for an
% omnidirectional mobile manipulator robot is analyzed in order to perform
% a pick and place task. We want the robot to move near the table without
% any collision, pick the object carefully and place it in the top of the
% moving platform. The desired motion is needed to be 10 seconds.

% In this file the trajectory planning is done in the JOINT SPACE
%% Clearing
clc;
clear;
close all;
tic;            % Start Clock to evaluate code
%% Initiallizing
t_start = 0;
ts=0.01;                    % Time step
t_end1 = 3;
t_plat=t_start:ts:t_end1;             % Time vector for platform motion

q0 = [2.6180 -0.6695 1.2719 3.1416 1.2002 -0.9821];% Arm initial conditions / Joints
q(1,:) = q0;
x1(1) = 0; % initial center position (x axis) in meters (m)
y1(1) = 0; % initial center position (y axis) in meters (m)
phi(1)= 0; % initial orientation in radians
qm(1,:) = [x1(1) y1(1) phi(1)];

% Table dimensions in meters (m)
Lax = 0.35;
Lay = 0.35;
Laz = 0.6;

% Object dimensions in meters (m)
ha = 0.1;
r = 0.025;
clr = [1 0 1];
nSides = 100;

% Moving Platform dimensions in meters (m)
Lbx = 0.75;
Lby = 1;
Lbz = 0.5;

% Object Initial position
p0A = [1.5 1.5 0.6];

% Distance between robot arm and platform's coordinates
pMB = [0 0.35 0.5];
dx = pMB(1);
dy = pMB(2);
dz = pMB(3);

% We want the robot to move in front of the table in a safe distance and
% carefully pick the object from the table
d_safe = 0.1;                     % 10 cm from the table as a safe distance

% OMR's distance from object
Dx = p0A(1);
Dy = p0A(2) - d_safe - Lay/2 - Lby/2;

% Tmb
TMB = [eye(3) [dx;dy;dz];0 0 0 1];
% T0M initially
T0M = [rotz(phi(1)) [x1(1);y1(1);0];0 0 0 1];
% T0B initially
T0B = T0M*TMB;

% Robot object for the robot arm
lwr = lwr_create();
% lwr.base = T0B;

% End-effector initial position
g = lwr.fkine(q0);
g = T0B*g.T;
hx(1) = g(1,4);     % End-effector initial x
hy(1) = g(2,4);     % End-effector initial y
hz(1) = g(3,4);     % End-effector initial z

%% Trajectories: JOINT SPACE

%% Trajectory 1: Moving the platform but not joints of the arm

% For the first trajectory we want only the platform to move without
% changing the orientation of the end-effector. Thus, for
% qe = [x y phi q0(1:6)] we want the angles of joints to stay the same,
% as well as the phi angle of rotation of the platform, and therefore
% moving changing only x and y of the platform. We start from the world
% frame and move the platform Dx and Dy distances on x and y axis
% accordingly. The desired motion is to start moving the platform and park
% it right in front of the target object in a safe distance.

start_plat = [0 0 0 q0];
finish_plat = [Dx Dy 0 q0];

% Calculating the positions and velocities of the platform and joints for
% the desired trajectory:
[q_1,q_1dot,c1] = jtrajMod(start_plat,finish_plat,t_plat);

% Calculating the position and velocity of the end-effector
k = 1;
for i = 1:length(q_1)
    if k == length(q_1)+1
        k = length(q_1)
    end
    % Update homogeneous transformations
    T0M = rt2tr(rotz(q_1(k,3)),[q_1(k,1);q_1(k,2);0]);
    T0B = T0M*TMB;
    % Jacobians throughout the trajectory
    Jp = [eye(2) [-dx*cos(q_1(k,3))-dy*sin(q_1(k,3));dx*sin(q_1(k,3))-dy*cos(q_1(k,3))];zeros(3);0 0 1];
    Ja = lwr.jacob0(q_1(k,4:9));
    J(:,:,k) = [Jp Ja];
    qd = q_1dot(k,:)';
    % Estimation of Ve (Velocity of the end-effector)
    Ve(:,k) = J(:,:,k)*qd;
    
    % Estimating the end-effector positions
    g = lwr.fkine(q_1(k,4:9));
    X1(:,:,k) = T0B*g.T;
    hx(k) = X1(1,4,k);
    hy(k) = X1(2,4,k);
    hz(k) = X1(3,4,k);
   
    k = k+1;
end
k = k-1;        % Update k to k = length(q_1)
%% Trajectory 2: Moving the joints of the arm but not the platform
% Pick and place operation:
% The trajectory is going to be a sum of 6 different motions:
% a) Part 1: Moving the arm above the object in a safe distance (0.3 m)
% to avoid collision / Ready-to-pick position.
% b) Part 2: Move the arm down to pick the object.
% c) Part 3: Move the arm up again in the same safe distance (0.3 m)
% d) Part 4 (Optional): Make a rotational motion moving only the first
% joint to avoid passing through singularity points e.g. extending the arm
% to move it to the other side. Stop the rotation at the antipodal point
% of the previous position.
% e) Part 5: Move the arm to a ready-to-place position above the {F} frame
% ina safe distance (0.3 m)
% d) Part 6: Placing the object to the {F} frame.
%
% For each of the motions, the desired orientation and position of the end
% effector is expressed in homogeneous transforms, then solving the inverse
% kinematics problem to get the angle positions at the desired via points,
% then use this angle positions as via points at the trajectory planning.
% When solving the inverse kinematics it is possible that there are more
% than one solutions -> angle positions to achieve the desired end-effector
% position could have more than one different solutions. To ensure that the
% arm avoids any collision we give a first estimation of the joints
% positions at the desired point, so we make sure the approach to the
% object is in a safe way from above. To get the estimations of such joint
% positions, lwr.teach method has been used while applying trial and error
% tests.
% At the end of the script and the whole motion we check if the
% arm passed through singularity points (det(J) = 0), to justify Part 4 of
% the motion.

% Pick and Place operation calculations:

% Distance of target object from the arm after the first motion.
dy_safe = d_safe+Lby/2-pMB(2)+Lay/2;

% Part 1:
TBup = rt2tr(g.R,[0;dy_safe;0.4]);
q_pickup = lwr.ikine(TBup,'q0',[pi/4 deg2rad(-77) deg2rad(41) 0 deg2rad(-62.5) 0]);

% Part 2:
TBA = rt2tr(g.R,[0;dy_safe;0.1]);
q_pick = lwr.ikine(TBA,'q0',[pi/4 deg2rad(-77) deg2rad(41) 0 deg2rad(-62.5) 0]);

% Part 3 has the same desired orientation and position with Part 1.

% Part 4:
TBFup1 = rt2tr(rotx(-pi)*rotz(-pi/2),[0;-dy_safe;0.4]);
q_placeup1 = lwr.ikine(TBFup1,'q0',[-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);

% Part 5:
TBFup2 = rt2tr(rotx(-pi)*rotz(-pi/2),[0;-0.35;0.3]);
q_placeup2 = lwr.ikine(TBFup2,'q0',[-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);

% Part 6:
TBF = rt2tr(rotx(-pi)*rotz(-pi/2),[0;-0.35;0]);
q_place = lwr.ikine(TBF,'q0',[-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);


% Move the arm to a ready-to-pick position
start_arm = [Dx Dy 0 q0];
f1 = [Dx Dy 0 q_pickup];
t1 = ts:ts:1.5;
[q_21,q_21dot,c21] = jtrajMod(start_arm,f1,t1);
[J21,Ja21,qd21,Ve21,X21,hx21,hy21,hz21,k] = mytrajLoop(q_21,q_21dot,hx,hy,hz,k,dx,dy,TMB,lwr);

% Move the arm to pick the object
f2 = [Dx Dy 0 q_pick];
t2 = ts:ts:1;
[q_22,q_22dot,c22] = jtrajMod(f1,f2,t2);
[J22,Ja22,qd22,Ve22,X22,hx22,hy22,hz22,k] = mytrajLoop(q_22,q_22dot,hx21,hy21,hz21,k,dx,dy,TMB,lwr);

% Move the arm up again, lifting the object
f3 = f1;
t3 = ts:ts:1;
[q_23,q_23dot,c23] = jtrajMod(f2,f3,t3);
[J23,Ja23,qd23,Ve23,X23,hx23,hy23,hz23,k] = mytrajLoop(q_23,q_23dot,hx22,hy22,hz22,k,dx,dy,TMB,lwr);

% Rotating the arm to avoid collision
f4 = [Dx Dy 0 q_placeup1];
t4 = ts:ts:1;
[q_24,q_24dot,c24] = jtrajMod(f3,f4,t4);
[J24,Ja24,qd24,Ve24,X24,hx24,hy24,hz24,k] = mytrajLoop(q_24,q_24dot,hx23,hy23,hz23,k,dx,dy,TMB,lwr);

% Moving to a ready-to-place position
f5 = [Dx Dy 0 q_placeup2];
t5 = ts:ts:1;
[q_25,q_25dot,c25] = jtrajMod(f4,f5,t5);
[J25,Ja25,qd25,Ve25,X25,hx25,hy25,hz25,k] = mytrajLoop(q_25,q_25dot,hx24,hy24,hz24,k,dx,dy,TMB,lwr);

% Finally, place the object to frame {F}
finish_arm = [Dx Dy 0 q_place];
t6 = ts:ts:1.5;
[q_26,q_26dot,c26] = jtrajMod(f5,finish_arm,t6);
[J26,Ja26,qd26,Ve26,X26,hx26,hy26,hz26,k] = mytrajLoop(q_26,q_26dot,hx25,hy25,hz25,k,dx,dy,TMB,lwr);

%% The whole trajectory: Trajectory 1 + Trajectory 2
% Via points
via = [finish_plat;f1;f2;f3;f4;f5;finish_arm];

% We take the q_whole with mstraj as it is the same as taking all of the
% jtrajMod functions used till now. To prove that, uncomment and run the
% lines 367-372 where we plot the motion of the arm in the Trajectory 2.
[q_whole,t] = mstraj(via,[],[3,1.5,1,1,1,1,1.5],start_plat,ts,[0.2 0.2 0.2 0.2 0.2 0.2 0.2]);

% Trajectory of the end-effector
hX = [hx';hx21';hx22';hx23';hx24';hx25';hx26'];
hY = [hy';hy21';hy22';hy23';hy24';hy25';hy26'];
hZ = [hz';hz21';hz22';hz23';hz24';hz25';hz26'];

% Velocities of platform and joints
q_wholedot = [q_1dot;q_21dot;q_22dot;q_23dot;q_24dot;q_25dot;q_26dot];
    
% End-effector velocities
Ve_whole = [Ve Ve21 Ve22 Ve23 Ve24 Ve25 Ve26];

% The whole Cartesian trajectory
X = cat(3,X1,X21,X22,X23,X24,X25,X26);

% Orientation of the end effector as Quaternion
Q = UnitQuaternion(X);

% Jacobian of the arm during the whole motion 
J_wholearm = cat(3,Ja,Ja21,Ja22,Ja23,Ja24,Ja25,Ja26);


toc;     % Stop Clock
%% Loop simulation
scene=figure;  % new figure
tam=get(0,'ScreenSize');
set(scene,'position',[tam(1) tam(2) tam(3) tam(4)]); % position and size figure in the screen
axis equal; % Set axis aspect ratios
axis([-3 3 -3 3 -0.5 1.5]); % Set axis limits 
grid on
% Draw the target cylinder-object.
target = plotCylinderWithCaps(r,p0A,ha,nSides,clr);
hold on
% Draw the initial position of the platform and the table
Plat = DrawCuboid([Lbx;Lby;Lbz],[0;0;Lbz/2]);
Table = DrawCuboid([Lax;Lay;Laz],[p0A(1);p0A(2);Laz/2]);

% Draw starting point of the end-effector
Point = plot3(hx(1),hy(1),hz(1),'rO','LineWidth',2); % Starting point
Trajectory1 = plot3(hx(1),hy(1),hz(1),'b','LineWidth',2); % Starting point

% Reset/Reinitialize the lwr.base -> Going to be updated to its true
% position inside the following loop.
lwr.base = eye(4);

% Orientation andposition of object/frame {A} of the object expressed in 
% the world frame.
T0A = transl(p0A);
fa = trplot(T0A,'frame','A','color',[1 0 1],'length',0.25);

% Draw world frame {0}
trplot(eye(4),'frame','0','color',[1 0 0],'length',0.25);

% Draw {B},{M},{F} frame in their initial positions
fb = trplot(lwr.base,'frame','B','color',[0 1 1],'length',0.25);
fm = trplot(eye(4),'frame','M','color',[0 1 1],'length',0.25);
ff = trplot(eye(4),'frame','F','color','k','length',0.25);

% Simulation loop
for k = 1:10:length(t)
    % Delete previous drawings to update them
    delete(Plat);
    delete(fb);
    delete(fm);
    delete(ff);
    axis([-3 3 -3 3 -0.5 1.5]); % Reset axis limits 
    
    % Update positions of platform
    x = q_whole(k,1);
    y = q_whole(k,2);
    
    % Update homogeneous transforms for platform frame {M}, base of
    % the arm frame {B} and frame of placing {F} related to world frame {0}
    g0b = [rotz(q_whole(k,3)) [x+pMB(1);y+pMB(2);pMB(3)];0 0 0 1];
    g0m = [rotz(q_whole(k,3)) [x;y;0];0 0 0 1];
    g0f = [rotz(q_whole(k,3)) [x;y;pMB(3)];0 0 0 1];
    
    % Update base of the lwr object
    lwr.base = g0b;
    
    % Show frames
    fb = trplot(g0b,'frame','B','color',[0 1 1],'length',0.25);
    fm = trplot(g0m,'frame','M','color',[1 1 0],'length',0.25);
    ff = trplot(g0f,'frame','F','color','k','length',0.25);
    
    % Show message if the robot finishes the platform-only move.
    if k > 3/ts
        text(-2,y,1.5,'- Robot Parked')
    end
    
    % Show message if the robot finishes the picking task.
    % Also, delete the object and frame {A} to symbolically show that the
    % object has been picked.
    if k >= 5.5/ts
        delete(target)
        delete(fa)
        text(-2,y,1.2,'- Object picked')
    end
    
    % Show message if the robot finishes the placing task.
    if k == length(t)-9
        text(-2,y,0.9,'- Object placed')
    end
    
    % Update platform and arm (using noname and nobase for the arm to fit
    % the mobile manipulator simulation)
    Plat = DrawCuboid([Lbx;Lby;Lbz],[x;y;Lbz/2]);
    lwr.plot(q_whole(k,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
    hold on
    % Update trajectory plot
    Trajectory1 = plot3(hX(1:k),hY(1:k),hZ(1:k),'b','LineWidth',2); % Starting point
    pause(ts)
    
    % Making it a gif
%     frame = getframe(scene); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if k == 1 
%         imwrite(imind,cm,'PickAndPlace.gif','gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,'PickAndPlace.gif','gif','WriteMode','append'); 
%     end 
end

% Re-plot the object in the position that is placed
target = plotCylinderWithCaps(r,[g0m(1,4);g0m(2,4);g0m(3,4)+Lbz],ha,nSides,clr);

% Pause to show that the motion has ended.
pause(3);
% Plot arm to the initial q0 position to show that the object has been
% placed and it is now on frame {f}
lwr.plot(q0,'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase');

% lwr.plot(q_21(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)
% lwr.plot(q_22(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)
% lwr.plot(q_23(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)
% lwr.plot(q_24(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)
% lwr.plot(q_25(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)
% lwr.plot(q_26(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','fps',20)

% Ensure that there were no singularity points at the whole trajectory:
for i = 1:size(J_wholearm,3)
    if det(J_wholearm(:,:,i)) == 0
        fprintf(['Singularity point at ' num2str(i) ' iteration.'])
    end
end

%% Plots of results

%% 1.Angles of joints / Platform and end-effector position -vs- time
% Angles of joints vs time
figure;
plot(t,q_whole(:,4:9))
title('Angles of joints -vs- Time')
xlabel('Time (s)')
ylabel('Angles (rad)')
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6')
grid on

% Platform position -vs- Time expressed in world frame
figure;
subplot(2,1,1)
plot(t,q_whole(:,1:3))
title('Platform coordinates (in world frame) -vs- Time')
xlabel('Time (s)')
ylabel('Position of the platform (m)')
legend('x_p_l_a_t_f_o_r_m','y_p_l_a_t_f_o_r_m','z_p_l_a_t_f_o_r_m')
grid on

% End-effector position -vs- Time expressed in world frame
tt = (0:ts:10)';
subplot(2,1,2)
plot(tt,hX,tt,hY,tt,hZ);
title('End-effector coordinates (in world frame) -vs- Time')
xlabel('Time (s)')
ylabel('Position of the end-effector (m)')
legend('x_e_e','y_e_e','z_e_e')
grid on

%% 2. 3D Representation of trajectories ( Platform / End-effector )
% Trajectory of end-effector and platform in 3D
sc = figure;
tam=get(0,'ScreenSize');
set(sc,'position',[tam(1) tam(2) tam(3) tam(4)]); % position and size figure in the screen
axis equal; % Set axis aspect ratios
axis([-3 3 -3 3 -0.5 1.5]); % Set axis limits 
grid on
plot3(hX,hY,hZ);
hold on
plot3(q_whole(:,1),q_whole(:,2),q_whole(:,3))
plot3(hX(1),hY(1),hZ(1),'rO')
plot3(q_whole(1,1),q_whole(1,2),q_whole(1,3),'kO')
title('Trajectory of the platform and the end-effector')
legend('Trajectory of the end-effector','Trajectory of the platform',...
    'Staring point of the end-effector','Starting point of the the platform')
grid on

%% 3. Angle of orientation of the end-effector as quaternion -vs- Time
% End-effector's orientation as Unit Quaternion
Q = UnitQuaternion(X);
[yaw,pitch,roll] = quat2angle(Q.double);
figure;
subplot(2,1,1)
plot(tt,yaw,tt,pitch,tt,roll)
title('RPY Angles of orientation of the end-effector -vs- Time')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Yaw','Pitch','Roll')
grid on

subplot(2,1,1)
Qd = Q.double;
thet = 2*atan2(1,Qd(:,1));
subplot(2,1,2)
plot(tt,thet)
title('Angle of unit quaternion of the orientation of the end effector -vs- Time')
xlabel('Time (s)')
ylabel('Angle (rad)')
grid on

%% 4. Velocities of platform / joints / end-effector -vs- Time
% Velocities of platform and joints
figure;
subplot(2,1,1)
plot(tt,q_wholedot(:,4:9))
title('Velocities of joints -vs- Time')
xlabel('Time (s)')
ylabel('Velocity (rad/s)')
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6')
grid on

subplot(2,1,2)
plot(tt,q_wholedot(:,1:3))
title('Platform velocity -vs- Time')
xlabel('Time (s)')
ylabel('Velocity of the platform (m)')
legend('ux_p_l_a_t_f_o_r_m','uy_p_l_a_t_f_o_r_m','uphi_p_l_a_t_f_o_r_m')
grid on

% Velocity of the end-effector
figure;
plot(tt,Ve_whole)
title('Velocity of the end-effector -vs- Time')
xlabel('Time (s)')
ylabel('End-effector velocity (m/s)')
legend('ux_e_e','uy_e_e','uz_e_e','wx_e_e','wy_e_e','wz_e_e')
grid on

%% 5. Coefficients of trajectories
C = [c1;c21;c22;c23;c24;c25;c26]
%% ------------------------- End of Project ------------------------------