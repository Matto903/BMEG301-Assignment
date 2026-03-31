clc, clear
% hello there my name is matthew

%Test From Oli 

gender = 2; % 1 for female, 2 for male
weight = 85; % Total body mass of subject in 'kg'
height = 184; % Total body height of subject in 'cm'

% Segment factors based of De lava Table
% Mass factor   Length factor   COM factor  Radius of Gyration factor

% Female factors
segment_table(1).factors = [ 
0.1478  368.5/1735  0.3612  0.369;  % Thigh
0.0481  432.3/1735  0.4416  0.291   % Shank

];
% Male factors
segment_table(2).factors = [
0.1416  422.2/1741  0.4095  0.329;  % Thigh
0.0433  434.0/1741  0.4459  0.255   % Shank
]; 

model(1).name = 'thigh'; model(1).color = 'k';
model(2).name = 'shank'; model(2).color = 'r';

table = sprintf('part Name\t Mass\t Length \t COM\t RGyration\t Inertia');
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
    
    t = sprintf('%s\t %.3f\t %.3f \t %.3f\t %.3f\t %.5f',model(part).name,...
        model(part).mass, model(part).length, model(part).com_local,...
        model(part).rgyration, model(part).inertia);
    disp(t)
end