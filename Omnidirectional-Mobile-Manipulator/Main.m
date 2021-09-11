%% ROBOTICS - Aristotle University of Thessaloniki 
% Project - Summer Semester 2020/2021
% Kavelidis Frantzis Dimitrios - AEM 9351 - kavelids@ece.auth.gr - ECE AUTH

%% Project Description:
% In this project the kinematics and trajectory planning for an
% omnidirectional mobile manipulator robot is analyzed in order to perform
% a pick and place task. We want the robot to move near the table without
% any collision, pick the object carefully and place it in the top of the
% moving platform. The desired motion is needed to be 10 seconds.

% In this file the trajectory planning is done in the TASK SPACE
%% Clearing
clc;
clear;
close all;
tic;        % Clock for code evaluation
%% Initiallizing
ts = 0.01;
tend = 20;

% Times for each segment
t_1 = (0:ts:5)';
t_21 = (0:ts:2.5)';
t_22 = (0:ts:2.5)';

t_23 = (0:ts:2.5)';
t_24 = (0:ts:3)';
t_25 = (0:ts:3)';
t_26 = (0:ts:1.5)';


q0 = [2.6180 -0.6695 1.2719 3.1416 1.2002 -0.9821];% Arm initial conditions / Joints


q(1,:) = q0;
x1(1) = 0; % initial center position (x axis) in meters (m)
y1(1) = 0; % initial center position (y axis) in meters (m)
phi(1)= 0; % initial orientation in radians
qm(1,:) = [x1(1) y1(1) phi(1)];

qe(1,:) = [qm q];
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

x(:,1) = zeros(6,1);

% Tmb
TMB = [eye(3) [dx;dy;dz];0 0 0 1];
% T0M initially
T0M = [rotz(phi(1)) [x1(1);y1(1);0];0 0 0 1];
% T0B initially
T0B = T0M*TMB;

% Robot object for the robot arm
lwr = lwr_create();
lwr.base = T0B;

% End-effector initial position
g = lwr.fkine(q0);
g1 = g.T;
hx(1) = g1(1,4);     % End-effector initial x
hy(1) = g1(2,4);     % End-effector initial y
hz(1) = g1(3,4);     % End-effector initial z
p(:,1) = [hx(1);hy(1);hz(1)];     
R(:,:,1) = g.R;                   % Initial end-effector orientation

%% Trajectories: TASK SPACE

%% Trajectory 1: Moving the platform but not joints of the arm

% For the first trajectory we want only the platform to move without
% changing the orientation of the end-effector. Thus we want the angles
% of joints to stay the same,
% as well as the phi angle of rotation of the platform, and therefore
% moving changing only x and y of end-effector/platform. We start from
% the world frame and move the platform Dx and Dy distances on x and y axis
% accordingly. The desired motion is to start moving the platform and park
% it right in front of the target object in a safe distance.

% End-effector trajectory planning:
Tstart1 = g.T;
[R1_start,p1_start] = tr2rt(Tstart1);% Initial end-effector pos and orientation
Q1_start = UnitQuaternion(R1_start);

Tfinish1 = rt2tr(g.R,[hx(1)+Dx;hy(1)+Dy;hz(1)]);
[R1_finish,p1_finish] = tr2rt(Tfinish1);
Q1_finish = UnitQuaternion(R1_finish);

% t_1 = (0:ts:3)';


% Zero velocities at the start and the end of the trajectory so we can park
% the robot first in front of the target.
qd1i = 0;
qd1f = 0;
[s1,sd1,sdd1,c1,cd1,cdd1] = tpolyMod(0,1,t_1,qd1i,qd1f);

% Interpolating poses by interpolating positions and orientations of the
% end-effector.
% Creating the trajectory p1(s) = p1_0 + s*(p1_1-p1_0)
dp1 = p1_finish-p1_start;
for i = 1:length(s1)
    p1(:,i) = p1_start + s1(i)*dp1;
    if i < length(s1)
        % Interpolating quaternions / orientations of the end-effector
        Q1(:,i) = Q1_start.interp(Q1_finish,s1(i));
    else
        Q1(:,i) = Q1_start.interp(Q1_finish,1);
    end
        % Velocity trajectory with interpolation
    p1d(:,i) = (dp1/norm(dp1))*sd1(i);
end
p1 = p1';
p1d = p1d';

% Quaternion to rpy
rpy1 = Q1.torpy;
rpy1d(1,:) = zeros(1,3);
for i = 1:(length(s1)-1)
    rpy1d(i+1,:) = (rpy1(i+1,:)-rpy1(i,:))/ts;
end

% Rebuild the homogeneous transformation matrices / get all the desired
% poses for the end-effector throughout the trajectory.
T_1 = rt2tr(rpy2r(rpy1(:,:)),p1(:,:)');

% End-effector's velocity.
Ve1 = [p1d rpy1d]';


% Our input is the vref = pinv(J)*Ve1 where J = [Jp Ja] and vref = [up u] 
% Since we want only the platform to move and we have a trajectory for the
% velocity of the end effector, we choose qdot to be 0 and qmdot to be same as
% the end-effector's velocity.

% Thus, the platform coordinates will be:
x1 = p1(:,1) - hx(1);
y1 = p1(:,2) - hy(1);

% First loop/ Trajectory 1 / Euler method 
for k = 1:length(s1)
    % Update homogeneous transformations
    T0M = rt2tr(rotz(phi(k)),[x1(k);y1(k);0]);
    T0B = T0M*TMB;
    lwr.base = T0B;
%     q(k,:) = lwr.ikine(T_1(:,:,k),'q0',qk(k,4:9));

    % Jacobians throughout the trajectory
    Jp = [eye(2) [-dx*cos(phi(k))-dy*sin(phi(k));dx*sin(phi(k))-dy*cos(phi(k))];zeros(3);0 0 1];
    Ja(:,:,k) = lwr.jacob0(q(k,:));
    J(:,:,k) = [Jp Ja(:,:,k)];
    
%     qe(k+1,:) = qe(k,:) + (pinv(J(:,:,k))*Ve1(:,k))';
%     qk(k+1,:) = qk(k,:) + (pinv(J(:,:,k))*(tr2delta(T_1(:,:,k))-tr2delta(lwr.fkine(qk(k,4:9)))))';
    pJ = pinv(J(:,:,k));
%     vref(:,k) = pJ*Ve1(:,k);

%   Our choice for the desired motion:
    vref(1:3,k) = pinv(Jp)*Ve1(:,k);

    % qmdot
    uf(k) = vref(1,k);        
    ul(k) = vref(2,k);        
    w(k) = vref(3,k);         
    
    qdot(k,:) = zeros(1,6);
    
    phi(k+1) = phi(k)+ts*w(k);
    q(k+1,:) = q(k,:)+ts*qdot(k,:); % Angular position
    
end
phi = phi(1:end-1)';
q = q(1:end-1,:);

% We get phi and w with almost zero -> 1e-13 significance,thus we assume they are 0;
phi = zeros(length(phi),1);
w = zeros(1,length(w));


fprintf('Trajectory 1: Calculated. \n')

%% Trajectory 2: Moving the joints of the arm but not the platform
% Pick and place operation:
% The trajectory is going to be a sum of 6 different motions:
% a) Part 1: Moving the arm above the object in a safe distance (0.3 m)
% to avoid collision / Ready-to-pick position.
% b) Part 2: Move the arm down to pick the object.
% c) Part 3: Move the arm up again in the same safe distance (0.3 m)
% d) Part 4: Move to a mid-point at the same height to avoid
% collision.
% e) Part 5: Move the arm to a ready-to-place position above the {F} frame
% in a safe distance (0.3 m)
% d) Part 6: Placing the object to the {F} frame.
%
% For each of the motions, the desired orientation and position of the end
% effector is expressed in homogeneous transforms, then interpolating poses
% between via points and solve the inverse kinematics at every point to get,
% the positions/angle of joints.
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

% % Distance of target object from the arm after the first motion.
dy_safe = d_safe+Lby/2-pMB(2)+Lay/2;
% % 
% % Part 1:
% TBup = rt2tr(roty(-pi)*rotz(-pi/4),[p0A(1);p0A(2);p0A(3)+0.3]);
% [R21_up,p21_up] = tr2rt(TBup);
% Q21_up = UnitQuaternion(R21_up);
Q21_up = Q1_finish;
p21_up = [p0A(1);p0A(2);p0A(3)+0.3];
TBup = rt2tr(R1_finish,p21_up);

% t_21 = (0:ts:1.5)';

% Zero velocities at the start and the end positions as we want to bring
% the arm in a ready-to-pick position above the target.
qd_21i = 0;
qd_21f = 0;
[T_21,Ve21,p21,pd21,Q21,rpy21,rpyd21,s21,sd21,sdd21,c21,cd21,cdd21] = ...
    myTraj(p1_finish,p21_up,Q1_finish,Q21_up,t_21,ts,qd_21i,qd_21f);

% Now the platform does not move, and thus we can consider that the
% manipulator now obeys the usual laws. Therefore:
for i = 1:length(s21)
    q(k,:) = lwr.ikine(T_21(:,:,i),'q0',...
    [pi/4 deg2rad(-77) deg2rad(41) 0 deg2rad(-62.5) 0]);    
    Ja21(:,:,i) = lwr.jacob0(q(k,:));
    vref21 = inv(Ja21(:,:,i))*Ve21(:,i);
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref21';
    k = k+1;
end
k = k-1;

fprintf('Trajectory 2.1: Calculated. \n')


% Part 2:
Q22_pick = Q1_finish;
p22_pick = p0A';
TBpick = rt2tr(Q22_pick.R,p22_pick);


% Zero velocities to stop motion at picking
qd22i = 0;
qd22f = 0;
% t_22 = (0:ts:1)';
[T_22,Ve22,p22,pd22,Q22,rpy22,rpyd22,s22,sd22,sdd22,c22,cd22,cdd22] = ...
    myTraj(p21_up,p22_pick,Q21_up,Q22_pick,t_22,ts,qd22i,qd22f);


for i = 1:length(s22)
    q(k,:) = lwr.ikine(T_22(:,:,i),'q0',...
    [pi/4 deg2rad(-77) deg2rad(41) 0 deg2rad(-62.5) 0]);    
    Ja22(:,:,i) = lwr.jacob0(q(k,:));
    vref22 = inv(Ja22(:,:,i))*Ve22(:,i);
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref22';
    k = k+1;
end
k = k-1;
fprintf('Trajectory 2.2: Calculated. \n')


% Part 3:
% t_23 = (0:ts:1)';

% Again to avoid any collision, we plan a discrete upwards motion after
% picking, so again choosing zero velocities:
qd23i = 0;
qd23f = 0;
[T_23,Ve23,p23,pd23,Q23,rpy23,rpyd23,s23,sd23,sdd23,c23,cd23,cdd23] = ...
    myTraj(p22_pick,p21_up,Q22_pick,Q21_up,t_23,ts,qd23i,qd23f);

for i = 1:length(s23)
    q(k,:) = lwr.ikine(T_23(:,:,i),'q0',...
    [pi/4 deg2rad(-77) deg2rad(41) 0 deg2rad(-62.5) 0]);    
    Ja23(:,:,i) = lwr.jacob0(q(k,:));
    vref23 = inv(Ja23(:,:,i))*Ve23(:,i);
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref23';
    k = k+1;
end
k = k-1;

fprintf('Trajectory 2.3: Calculated. \n')


% Part 4:
g = lwr.base.T;
gx = g(1,4);
gy = g(2,4);
TBFmid = rt2tr(rotx(-pi),[gx+dy_safe;gy;p0A(3)+0.3]);
[R24_mid,p24_mid] = tr2rt(TBFmid);
Q24_mid = UnitQuaternion(R24_mid);

qd24i = 0;
qd24f = 0;
% t_24 = (0:ts:1)';
[T_24,Ve24,p24,pd24,Q24,rpy24,rpyd24,s24,sd24,sdd24,c24,cd24,cdd24] = ...
    myTraj(p21_up,p24_mid,Q21_up,Q24_mid,t_24,ts,qd24i,qd24f);

% Fixing
% Ve24(4,250) = (Ve24(4,251)+Ve24(4,249))/2;
% Ve24(6,83) = (Ve24(6,82)+Ve24(6,84))/2;


for i = 1:length(s24)
    q(k,:) = lwr.ikine(T_24(:,:,i),'q0',...
    [-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);    
    Ja24(:,:,i) = lwr.jacob0(q(k,:));
    vref24 = inv(Ja24(:,:,i))*Ve24(:,i);
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref24';
    k = k+1;
end
k = k-1;

via = [gx+dy_safe gy p0A(3)+0.3;gx gy-pMB(2) pMB(3)+0.3];
start45 = [p0A(1);p0A(2);p0A(3)+0.3]';
p45 = mstraj(via,[],[1 1],start45,ts,.75);

fprintf('Trajectory 2.4: Calculated. \n')


% Part 5:
TBFup = rt2tr(rotx(-pi)*rotz(-pi/2),[gx;gy-pMB(2);pMB(3)+0.3]);
[R25_up,p25_up] = tr2rt(TBFup);
Q25_up = UnitQuaternion(R25_up);

qd25i = 0;
qd25f = 0;
% t_25 = (0:ts:1)';
[T_25,Ve25,p25,pd25,Q25,rpy25,rpyd25,s25,sd25,sdd25,c25,cd25,cdd25] = ...
    myTraj(p24_mid,p25_up,Q24_mid,Q25_up,t_25,ts,qd25i,qd25f);

% qk = zeros(length(s25),6);
% qk(1,:) = [-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0];
for i = 1:length(s25)
    q(k,:) = lwr.ikine(T_25(:,:,i),'q0',[-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);
    Ja25(:,:,i) = lwr.jacob0(q(k,:));
    vref25 = inv(Ja25(:,:,i))*Ve25(:,i);
%     qk(i+1,:) = [-(i/length(s25))*pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref25';
    k = k+1;
end
k = k-1;

fprintf('Trajectory 2.5: Calculated. \n')

% Part 6:

TBFplace = rt2tr(rotx(-pi)*rotz(-pi/2),[gx;T0M(2,4);pMB(3)]);
[R26_place,p26_place] = tr2rt(TBFplace);
Q26_place = UnitQuaternion(R26_place);

qd26i = 0;
qd26f = 0;
% t_26 = (0:ts:1.5)';
[T_26,Ve26,p26,pd26,Q26,rpy26,rpyd26,s26,sd26,sdd26,c26,cd26,cdd26] = ...
    myTraj(p25_up,p26_place,Q25_up,Q26_place,t_26,ts,qd26i,qd26f);

for i = 1:length(s26)
    q(k,:) = lwr.ikine(T_26(:,:,i),'q0',[-pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]);
    Ja26(:,:,i) = lwr.jacob0(q(k,:));
    vref26 = inv(Ja26(:,:,i))*Ve26(:,i);
%     qk(i+1,:) = [-(i/length(s25))*pi/2 deg2rad(-53) deg2rad(125) 0 deg2rad(-2.5) 0]
    uf(k) = 0;
    ul(k) = 0;
    w(k) = 0;
    x1(k) = x1(length(p1)-1);
    y1(k) = y1(length(p1)-1);
    phi(k) = phi(length(p1)-1);
    qdot(k,:) = vref26';
    k = k+1;
end
k = k-1;

fprintf('Trajectory 2.6: Calculated. \n')

%% Whole trajectory
% via = [p1_finish';p21_up';p22_pick';p21_up';p24_mid';p25_up';p26_place']
% p = mstraj(via,[],[5,2.5,2.5,2.5,3,3,1.5],p1_start',ts,0.05)

% Joints
q_whole = [x1 y1 phi q];

% Cartesian Coordinates of the end-effector
p = [p1(1:end-1,:);p21(1:end-1,:);p22(1:end-1,:);p23(1:end-1,:);p24(1:end-1,:);p25(1:end-1,:);p26];
% p = [p1(1:end-1,:);p21(1:end-1,:);p22(1:end-1,:);p23(1:end-1,:);p45;p26];

% Whole cartesian trajectory as homogeneous transform
T = cat(3,T_1(:,:,1:end-1),T_21(:,:,1:end-1),T_22(:,:,1:end-1),...
    T_23(:,:,1:end-1),T_24(:,:,1:end-1),T_25(:,:,1:end-1),T_26);

% Quaternion representation of end-effector's orientation
Q = UnitQuaternion(T);

% Jacobian of the arm in the whole movement
Ja_whole = cat(3,Ja(:,:,1:end-1),Ja21(:,:,1:end-1),Ja22(:,:,1:end-1),...
    Ja23(:,:,1:end-1),Ja24(:,:,1:end-1),Ja25(:,:,1:end-1),Ja26);

% Velocities:
qd_whole = [uf' ul' w' qdot];

% Velocity of the end effector
Ve_whole = [Ve1(:,1:end-1)';Ve21(:,1:end-1)';Ve22(:,1:end-1)';...
    Ve23(:,1:end-1)';Ve24(:,1:end-1)';Ve25(:,1:end-1)';Ve26'];

toc;        % Stop clock
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
t =(0:ts:tend)';
for k = 1:20:length(t)
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
    if k > 5/ts
        text(-2,y,1.5,'- Robot Parked')
    end
    
    % Show message if the robot finishes the picking task.
    % Also, delete the object and frame {A} to symbolically show that the
    % object has been picked.
    
    if k >= 10/ts
        delete(target)
        delete(fa)
        text(-2,y,1.2,'- Object picked')
    end
%     
    % Show message if the robot finishes the placing task.
    if k == length(t)
        text(-2,y,0.9,'- Object placed')
    end
    
    % Update platform and arm (using noname and nobase for the arm to fit
    % the mobile manipulator simulation)
    Plat = DrawCuboid([Lbx;Lby;Lbz],[x;y;Lbz/2]);
    lwr.plot(q(k,:),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','floorlevel',0)
    
    hold on
    
    % Update trajectory plot
    Trajectory1 = plot3(p(1:k,1),p(1:k,2),p(1:k,3),'b','LineWidth',2); % Starting point
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
% 
% Re-plot the object in the position that is placed
target = plotCylinderWithCaps(r,[g0m(1,4);g0m(2,4);g0m(3,4)+Lbz],ha,nSides,clr);

% Pause to show that the motion has ended.
pause(3);
% Plot arm to the initial q0 position to show that the object has been
% placed and it is now on frame {f}
lwr.plot(q0,'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase','floorlevel',0);

% % lwr.plot(q_21(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% % lwr.plot(q_22(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% % delete(target);
% % lwr.plot(q_23(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% % lwr.plot(q_24(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% % lwr.plot(q_25(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% % lwr.plot(q_26(:,4:9),'noname','workspace',[-3 3 -3 3 -0.5 1.5],'nobase')
% 
% % Ensure that there were no singularity points at the whole trajectory:
% for i = 1:size(J_wholearm,3)
%     if det(J_wholearm(:,:,i)) == 0
%         fprintf(['Singularity point at ' num2str(i) ' iteration.'])
%     end
% end
% 
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
subplot(2,1,2)
plot(t,p);
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
plot3(p(:,1),p(:,2),p(:,3))
hold on
plot3(q_whole(:,1),q_whole(:,2),q_whole(:,3))
plot3(p(1,1),p(1,2),p(1,3),'rO')
plot3(q_whole(1,1),q_whole(1,2),q_whole(1,3),'kO')
title('Trajectory of the platform and the end-effector')
legend('Trajectory of the end-effector','Trajectory of the platform',...
    'Staring point of the end-effector','Starting point of the the platform')
grid on

%% 3. Angle of orientation of the end-effector as quaternion -vs- Time
% End-effector's orientation as Unit Quaternion
[yaw,pitch,roll] = quat2angle(Q.double);
figure;
subplot(2,1,1)
plot(t,yaw,t,pitch,t,roll)
title('RPY Angles of orientation of the end-effector -vs- Time')
xlabel('Time (s)')
ylabel('Angle (rad)')
legend('Yaw','Pitch','Roll')
grid on

subplot(2,1,1)
Qd = Q.double;
thet = 2*atan2(1,Qd(:,1));
subplot(2,1,2)
plot(t,thet)
title('Angle of unit quaternion of the orientation of the end effector -vs- Time')
xlabel('Time (s)')
ylabel('Angle (rad)')
grid on

%% 4. Velocities of platform / joints / end-effector -vs- Time
% Velocities of platform and joints
figure;
subplot(2,1,1)
plot(t,qd_whole(:,4:9))
title('Velocities of joints -vs- Time')
xlabel('Time (s)')
ylabel('Velocity (rad/s)')
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6')
% axis([0 20 -20 10])
grid on

subplot(2,1,2)
plot(t,qd_whole(:,1:3))
title('Platform velocity -vs- Time')
xlabel('Time (s)')
ylabel('Velocity of the platform (m)')
legend('ux_p_l_a_t_f_o_r_m','uy_p_l_a_t_f_o_r_m','uphi_p_l_a_t_f_o_r_m')
% axis([0 10 -10 20])
grid on

% Velocity of the end-effector
figure;
plot(t,Ve_whole)
title('Velocity of the end-effector -vs- Time')
xlabel('Time (s)')
ylabel('End-effector velocity (m/s)')
legend('ux_e_e','uy_e_e','uz_e_e','wx_e_e','wy_e_e','wz_e_e')
% axis([0 20 -10 20])
grid on

%% 5. Trajectory Coefficients
C = [c1;c21;c22;c23;c24;c25;c26];
C = flip(C,2)
Cd = [cd1;cd21;cd22;cd23;cd24;cd25;cd26];
Cd = flip(Cd,2)


%% ------------------------- End of Project ------------------------------