 %% BMEG301 Reverse Dynamic Solution For Device 3
% Authors: Bailey Phillips, Matthew Nuske, Oliver Michels

%% Assumptions

% The thigh is consider to be pivoted at the hip joint and has one DoF in flexion-extension. 
% Distal limb (shank) is fixed such that the longitudinal axis line up. 
% All other limbs can be ignored.

clc, clear

%% Parameters

% Subject Paramters
gender = 2; % 1 for female, 2 for male
weight = 85; % Total body mass of subject in 'kg'
height = 1.84; % Total body height of subject in 'm'


%% Segment Factor Tables

% Segment factors based of De lava Table
% [Mass factor   Length factor   COM factor  Radius of Gyration factor]

% Female factors
segment_table(1).factors = [ 
0.1478  368.5/1735  0.3612  0.369;  % Thigh (female)
0.0481  432.3/1735  0.4416  0.271   % Shank (female)
];

% Male factors
segment_table(2).factors = [
0.1416  422.2/1741  0.4095  0.329;  % Thigh (male)
0.0433  434.0/1741  0.4459  0.255   % Shank (male)
];  


%% TASK 1 
% Create Structures With All of The Model's Segment Details 

disp('Human only system:');

model(1).name = 'Thigh'; model(1).color = 'k';
model(2).name = 'Shank'; model(2).color = 'r';

table = sprintf('Part Name \t Mass (kg)\t Length (m)\t COM (m)\t RGyration (m)\t Inertia (kg.m^2)');
disp([newline, table]);

% Assumption: The distal limb (shank) is fixed such that longitudinal axes
% line up, thus all other limbs can be ignored
for part = 1:2

    model(part).mass = weight*segment_table(gender).factors(part,1);
    model(part).length = height*segment_table(gender).factors(part,2);
    model(part).com_local = model(part).length*segment_table(gender).factors(part,3);

   if part < 2
        model(part).com_from_O = model(part).com_local;
    else
        model(part).com_from_O = model(part).com_local + sum([model(1:part-1).length]);
    end
    model(part).rgyration = model(part).length*segment_table(gender).factors(part,4);
    model(part).inertia = model(part).mass*(model(part).rgyration).^2;
    
    t = sprintf('%s\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.5f',model(part).name,...
        model(part).mass, model(part).length, model(part).com_local,...
        model(part).rgyration, model(part).inertia);
    disp(t);
end

% System Center of Mass
mass_system = sum([model(:).mass]);

% display system mass
disp([newline, 'System mass = ', num2str(mass_system), ' kg']);

r_num = [model(:).mass].*[model(:).com_from_O];
CoM_system = sum(r_num)/mass_system;
disp(['System center of mass position = ',num2str(CoM_system),' m']);

% System Inertia About The Hip
I_seg_o = [model(:).inertia] + [model(:).mass].*([model(:).com_from_O].^2);
Inertia_system = sum(I_seg_o);
disp(['System mass moment of inertia at O = ',num2str(Inertia_system),' kg.m^2']);

%% Simulation Parameters
fps = 50;                           % Frames per second 
                    % Motion duration in seconds

% These are the first things i saw for the rom of the thigh about the hip.
% we can change this if we find its different with more research - Oli
theta_min = -30 * (pi/180);           % Minimum angle at extension in rad
theta_max = 120 * (pi/180);         % Maximum angle at flexion in rad

%% TASK 2 & 3
durations=[5,3,1];

t_store = cell(length(durations),1);
theta_store = cell(length(durations),1);
omega_store = cell(length(durations),1);
alpha_store = cell(length(durations),1);
M_store = cell(length(durations),1);
P_store = cell(length(durations),1);

for i=1:length(durations)
    T_motion = durations(i);   
    % Simulation 
    Num_Frames = ceil(T_motion*fps);                        % Calculates number of frames in the simulation
    T_simulation = linspace(0, T_motion, Num_Frames);       % Creates a linearly spaced time vector for motion duration

    % Sigmoid Function To Calculate Joint Angle
    k = 10 / T_motion;                                                      % Steepness of sigmoid function (scales with the motion duration)
    s = 1 ./ (1 + exp(-k .* (T_simulation - (T_motion/2))));                % Sigmoid function
    theta = (theta_min + (theta_max - theta_min) .* s) .* (180/pi);         % Joint angle in degrees
    omega = gradient(theta,T_simulation);                                   % Angular velocity
    alpha = gradient(omega,T_simulation);                                   % Angular acceleration
    g = -9.81; % Gravity
    M_hip = (Inertia_system * alpha) - (g * CoM_system * mass_system * cosd(theta)); % Moment about Hip in Nm
    P = M_hip .* omega; % Joint power Watts

    t_store{i} = T_simulation;
    theta_store{i} = theta;
    omega_store{i} = omega;
    alpha_store{i} = alpha;
    M_store{i} = M_hip;
    P_store{i} = P;
end
%% Plotting Results
figure;
label = {'Slow ', 'Normal ', 'Fast '};

for i = 1:3
    % Joint Angle Over Time
    subplot(3, 3, i);
    plot(t_store{i}, theta_store{i});
    ylabel('Angle (deg)');
    xlabel('Time (s)')
    title([label{i} 'Joint Angle']);
    
    % Joint Angular Velocity Over Time
    subplot(3, 3, i + 3);
    plot(t_store{i}, omega_store{i})
    ylabel('Velocity (deg/s)');
    xlabel('Time (s)')
    title([label{i} 'Joint Anglular Velocity']);
    
    % Joint Angular Accelleration Over Time
    subplot(3, 3, i + 6);
    plot(t_store{i}, alpha_store{i});
    ylabel('Accelleration (deg/s^2)');
    xlabel('Time (s)')
    title([label{i} 'Joint Anglular Accelleration']);
end

figure;
for i = 1:3
    % Joint Moment Over Time
    subplot(2, 3, i);
    plot(t_store{i}, M_store{i});
    ylabel('Angle (deg)');
    xlabel('Time (s)')
    title([label{i} 'Joint Moment']);
    
    % Joint Power Over Time
    subplot(2, 3, i + 3);
    plot(t_store{i}, P_store{i})
    ylabel('Power (W)');
    xlabel('Time (s)')
    title([label{i} 'Joint Power']);
    
end
hold off;


%% Part 5

mass_device_total = 8;         % EXAMPLE device mass in kg
mass_device_pivoted = 0.6 .* mass_device_total;

model_device = model;

model_device(1).rgyration = model(1).rgyration;
model_device(1).mass = model_device(1).mass + mass_device_pivoted;
model_device(1).inertia = model_device(1).mass * (model_device(1).rgyration)^2;

disp([newline, 'Combined System:']);
disp([newline, table]);

for part = 1:2
    t_device = sprintf('%s\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.5f',model_device(part).name,...
            model_device(part).mass, model_device(part).length, model_device(part).com_local,...
            model_device(part).rgyration, model_device(part).inertia);
        disp(t_device);
end


% Combined System Center of Mass
mass_system_device = sum([model_device(:).mass]);

% Display combined system mass
disp([newline, 'Combined system mass = ', num2str(mass_system_device), ' kg']);

r_num_device = [model_device(:).mass] .* [model_device(:).com_from_O];
CoM_system_device = sum(r_num_device)/mass_system_device;
disp(['Combined system center of mass position = ',num2str(CoM_system_device),' m']);

% Combined System Inertia About The Hip
I_seg_o_device = [model_device(:).inertia] + [model_device(:).mass].*([model_device(:).com_from_O].^2);
Inertia_system_device = sum(I_seg_o_device);
disp(['Combined system mass moment of inertia at O = ',num2str(Inertia_system_device),' kg.m^2']);


M_hip_system_device = (Inertia_system_device .* alpha_store{3}) - (g * CoM_system_device * mass_system_device * cosd(theta_store{3})); % Moment about Hip in Nm
P_system_device = M_hip_system_device .* omega_store{3}; % Joint power Watts

P_device = P_system_device - P_store{3};
P_required = max(P_device + (0.5 .* P_store{3}));






%% Animation
caseNum = 1;              % 1 = slow, 2 = normal, 3 = fast
t1 = t_store{caseNum};
theta1 = theta_store{caseNum};

L = model(1).length + model(2).length;

figure;

hLeg = plot([0 0], [0 -L], 'LineWidth', 4);
hold on
plot(0,0,'ko','MarkerFaceColor','k')
axis equal
axis([-L L -L L])
grid on
xlabel('X Position (m)')
ylabel('Y Position (m)')

for j = 1:length(t1)
    x = L * sind(theta1(j));
    y = -L * cosd(theta1(j));

    set(hLeg, 'XData', [0 x], 'YData', [0 y])
    
title(sprintf('Sagittal Plane Motion (%s), t = %.2f s', label{caseNum}, t1(j)))

    drawnow
    pause(0.02)
    hold off;
end
