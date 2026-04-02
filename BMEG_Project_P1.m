 %% BMEG301 Reverse Dynamic Solution For Device 3
% Authors: Bailey Phillips, Matthew Nuske, Oliver Michels

%% Assumptions

% The thigh is consider to be pivoted at the hip joint and has one DoF in flexion-extension. 
% Distal limb (shank) is fixed such that the longitudinal axis line up. 
% All other limbs can be ignored.
% Pivoted part of the device contains 60% of device mass
% Non pivoted part of the device is rigidly attatched to the human base limb (hip)
% Pivot limb and device are coincident
% Device does not change radius of gyration
% Reasonable values for moment arm and surface area of the device attatching strap

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

% Start displaying human system information
disp('Human only system:');
model(1).name = 'Thigh'; model(1).color = 'k';
model(2).name = 'Shank'; model(2).color = 'r';
table = sprintf('Part Name \t Mass (kg)\t Length (m)\t COM (m)\t RGyration (m)\t Inertia (kg.m^2)');
disp([newline, table]);


% Create structures with all of the model's segment details 
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
    
    % Finish displaying human system
    t = sprintf('%s\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.5f',model(part).name,...
        model(part).mass, model(part).length, model(part).com_local,...
        model(part).rgyration, model(part).inertia);
    disp(t);
end

mass_system = sum([model(:).mass]);         % System center of Mass

% display system mass
disp([newline, 'System mass = ', num2str(mass_system), ' kg']);

r_num = [model(:).mass].*[model(:).com_from_O];
CoM_system = sum(r_num)/mass_system;
disp(['System center of mass position = ',num2str(CoM_system),' m']);

% System inertia about the hip
I_seg_o = [model(:).inertia] + [model(:).mass].*([model(:).com_from_O].^2);
Inertia_system = sum(I_seg_o);
disp(['System mass moment of inertia at O = ',num2str(Inertia_system),' kg.m^2']);

%% Simulation Parameters

fps = 50;                       % Frames per second 
theta_min = -30 * (pi/180);     % Minimum angle at extension in rad
theta_max = 120 * (pi/180);     % Maximum angle at flexion in rad

%% TASK 2 & 3

durations=[5,3,2]; % Motion durations [slow, medium, fast]

% Create cells to store information for each motion duration separately
t_store = cell(length(durations),1);        % Time storage cell
theta_store = cell(length(durations),1);    % Theta storage cell
omega_store = cell(length(durations),1);    % Omega storage cell
alpha_store = cell(length(durations),1);    % Alpha storage cell
M_store = cell(length(durations),1);        % Moment storage cell
P_store = cell(length(durations),1);        % Power storage cell

for i=1:length(durations)
    T_motion = durations(i);   % Motion duration in seconds
 
    Num_Frames = ceil(T_motion*fps);                        % Calculates number of frames in the simulation
    T_simulation = linspace(0, T_motion, Num_Frames);       % Creates a linearly spaced time vector for motion duration

    % Create normalised sigmoid function to calculate joint angle across time
    k = 0.1 / T_motion;                                                 % Steepness of sigmoid function
    s_raw = 1 ./ (1 + exp(-k .* (T_simulation - (T_motion/2))));        % Raw function to normalise s between 0-1
    s = (s_raw - s_raw(1)) ./ (s_raw(end) - s_raw(1));                  % Normalised sigmoid function

    % Apply sigmoid function to determine joint angle across time, and calulate angular velocity and acceleration
    theta = (theta_min + (theta_max - theta_min) .* s);      % Joint angle in rad
    omega = gradient(theta,T_simulation);                   % Angular velocity rad/s
    alpha = gradient(omega,T_simulation);                   % Angular  acceleration rad/s^2
    
    % Calculate the moment and power about the hip
    g = 9.81;                                                                           % Gravity
    M_hip = (Inertia_system * alpha) + (g * CoM_system * mass_system * sin(theta));     % Moment about Hip in Nm
    P = M_hip .* omega;                                                                 % Joint power Watts

    % Store information for respective motion duration
    t_store{i} = T_simulation;
    theta_store{i} = theta;
    omega_store{i} = omega;
    alpha_store{i} = alpha;
    M_store{i} = M_hip;
    P_store{i} = P;
end

%% Task 2 Plotting 

figure;
sgtitle('Task 2: Human System');
label = {'Slow ', 'Normal ', 'Fast '};

for i = 1:3
    % Joint Angle Over Time
    subplot(3, 3, i);
    plot(t_store{i}, theta_store{i} .* (180/pi));
    ylabel('Angle (deg)');
    xlabel('Time (s)')
    title([label{i} 'Joint Angle']);
    
    % Joint Angular Velocity Over Time
    subplot(3, 3, i + 3);
    plot(t_store{i}, omega_store{i} .* (180/pi))
    ylabel('Velocity (deg/s)');
    xlabel('Time (s)')
    title([label{i} 'Joint Anglular Velocity']);
    
    % Joint Angular   Acceleration Over Time
    subplot(3, 3, i + 6);
    plot(t_store{i}, alpha_store{i} .* (180/pi));
    ylabel('  Acceleration (deg/s^2)');
    xlabel('Time (s)')
    title([label{i} 'Joint Anglular   Acceleration']);
end

%% Task 3 Plotting

figure;
sgtitle('Task 3: Human System');

for i = 1:3
    % Joint Moment Over Time
    subplot(2, 3, i);
    plot(t_store{i}, M_store{i});
    ylabel('Moment (Nm)');
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



%% Tasks 4 & 5

mass_device_total = 13;                             % Device mass in kg 
mass_device_pivoted = 0.6 .* mass_device_total;     % Pivoted device mass in kg

% Define a sctructure for the human-device model
model_device = model;

model_device(1).rgyration = model(1).rgyration;                                     % Radius of gyration remains the same (as stated in assumptions)
model_device(1).mass = model_device(1).mass + mass_device_pivoted;                  % Mass of the system is taken to be the original mass, increased by the pivoted mass of the device
model_device(1).inertia = model_device(1).mass * (model_device(1).rgyration)^2;     % Calculate the inertia of the human-device system

% Display new human-device system information
disp([newline, 'Human-Device System:']);
disp([newline, table]);

for part = 1:2
    t_device = sprintf('%s\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t\t %.5f',model_device(part).name,...
            model_device(part).mass, model_device(part).length, model_device(part).com_local,...
            model_device(part).rgyration, model_device(part).inertia);
        disp(t_device);
end


% Calculate the human-device system center of mass
mass_system_device = sum([model_device(:).mass]);

% Display human-device system mass
disp([newline, 'Human-Device system mass = ', num2str(mass_system_device), ' kg']);

r_num_device = [model_device(:).mass] .* [model_device(:).com_from_O];
CoM_system_device = sum(r_num_device)/mass_system_device;
disp(['Human-Device system center of mass position = ',num2str(CoM_system_device),' m']);

% Calculate and display the human-device system inertia about the hip
I_seg_o_device = [model_device(:).inertia] + [model_device(:).mass].*([model_device(:).com_from_O].^2);
Inertia_system_device = sum(I_seg_o_device);

disp(['Human-Device system mass moment of inertia at O = ',num2str(Inertia_system_device),' kg.m^2']);

% Create cells to store human-device moment and power for each motion duration separately
M_system_device_store = cell(length(durations),1);      % Human-device system moment store cell
P_system_device_store = cell(length(durations),1);      % Human-device system power store cell

for i=1:length(durations)
    % Calculate the moment and power of the system about the hip
    M_hip_system_device = (Inertia_system_device * alpha_store{i}) + (g * CoM_system_device * mass_system_device * sin(theta_store{i}));    % Moment about Hip in Nm
    P_system_device = M_hip_system_device .* omega_store{i};                                                                                % Joint power Watts

    % Store information for respective motion duration
    M_system_device_store{i} = M_hip_system_device;   
    P_system_device_store{i} = P_system_device;
end

% Calculate maximum power based upon provided formula and parameters
a = 10.0;
b = 50.0;
c = 1.0;

P_max = (exp(mass_device_total - a) ./ (1 + exp(mass_device_total - a))) * b + c * mass_device_total        % Maximum power from the device

% Determine the pressure on the thigh from the device
M_Device = M_hip_system_device - M_store{3};


% Assuming the thigh to be a perfect cylinder, calculate the thigh circumference
thigh_density = 1050;                               % Thigh density value (from research) in kg/m^3
thigh_volume = model(1).mass / thigh_density;       % Thigh volume in m^3
thigh_cs_area = thigh_volume / model(1).length;     % Thigh cross sectional area in m^2
thigh_radius = sqrt(thigh_cs_area / pi);            % Radius of thigh in m
thigh_circumference = 2 * pi * thigh_radius;        % Circumference of the thigh, cylinder in m

% Assumed device cuff geometry
r_arm = 0.3;                                        % Assumed effective moment arm from hip to thigh cuff in m
strap_width = 0.08;                                 % Assumed width of strap in m
contact_length = 0.5 * thigh_circumference;         % Contact length assuming half of thigh circumference is in contact in m
A_strap = strap_width * contact_length;             % Area of the strap in m^2

% Strap force magnitude
F_strap = abs(M_Device) ./ r_arm;         % N

% Pressure on thigh
pressure = F_strap ./ A_strap;            % Pa
pressure_kPa = pressure / 1000;           % kPa

% Peak values
F_peak = max(F_strap);
pressure_peak = max(pressure);
pressure_peak_kPa = max(pressure_kPa)
pressure_avg = mean(pressure);          % Pa
pressure_avg_kPa = mean(pressure_kPa)  % kPa


%% Task 4 Plotting

% Calulate the respective power for a range of device masses
m_device = linspace(0, 20, 500);                                                % 500 linear spaced values 0-20
P_max1 = (exp(m_device - a) ./ (1 + exp(m_device - a))) * b + c * m_device;     % 500 corrosponding power values

% Plot the power to mass ratio function
figure;
sgtitle('Task 4: Power to Mass Function')
plot(m_device, P_max1, 'LineWidth', 1.5)

xlabel('Total Device Mass (kg)')

grid on

%% Task 5 Plotting

figure;
sgtitle('Task 5: Human + Device System');

for i = 1:3
    % Joint Moment Over Time
    subplot(2, 3, i);
    plot(t_store{i}, M_system_device_store{i});
    ylabel('Moment (Nm)');
    xlabel('Time (s)')
    title([label{i} 'Joint Moment']);
    
    % Joint Power Over Time
    subplot(2, 3, i + 3);
    plot(t_store{i}, P_system_device_store{i})
    ylabel('Power (W)');
    xlabel('Time (s)')
    title([label{i} 'Joint Power']);
    
end

P_device = P_system_device - P_store{3};
P_required = max(P_device + (0.5 .* P_store{3}))

%% Animation
caseNum = 3;              % 1 = slow, 2 = normal, 3 = fast
t1 = t_store{caseNum};
theta1 = theta_store{caseNum};

L = model(1).length + model(2).length;

figure;
sgtitle('Task 3: Human Limb Animation');
hLeg = plot([0 0], [0 -L], 'LineWidth', 4);
hold on
plot(0,0,'ko','MarkerFaceColor','k')
axis equal
axis([-L L -L L])
grid on
xlabel('X Position (m)')
ylabel('Y Position (m)')

for j = 1:length(t1)
    x = L * sin(theta1(j));
    y = -L * cos(theta1(j));

    set(hLeg, 'XData', [0 x], 'YData', [0 y])
    
title(sprintf('Sagittal Plane Motion (%s), t = %.2f s', label{caseNum}, t1(j)))

    drawnow
    pause(0.02)
    hold off;
end
