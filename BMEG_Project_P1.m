clc, clear

gender = 2; % 1 for female, 2 for male
weight = 85; % Total body mass of subject in 'kg'
height = 1.84; % Total body height of subject in 'm'

% Segment factors based of De lava Table
% Mass factor   Length factor   COM factor  Radius of Gyration factor

% Female factors
segment_table(1).factors = [ 
0.1478  368.5/1735  0.3612  0.369;  % Thigh
0.0481  432.3/1735  0.4416  0.271   % Shank
];

% Male factors
segment_table(2).factors = [
0.1416  422.2/1741  0.4095  0.329;  % Thigh
0.0433  434.0/1741  0.4459  0.255   % Shank
]; 

model(1).name = 'Thigh'; model(1).color = 'k';
model(2).name = 'Shank'; model(2).color = 'r';

table = sprintf('Part Name\t Mass\t\t Length \t COM\t\t RGyration\t Inertia');
disp(table)

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
    disp(t)
end

% System center of mass
mass_system = sum([model(:).mass]);
r_num = [model(:).mass].*[model(:).com_from_O];
CoM_system = sum(r_num)/mass_system;
disp([newline, 'System center of mass position = ',num2str(CoM_system),' m']);

% System inertia about the hip
I_seg_o = [model(:).inertia] + [model(:).mass].*([model(:).com_from_O].^2);
Inertia_system = sum(I_seg_o);
disp(['System mass moment of inertia at O = ',num2str(Inertia_system),' kg.m^2']);
