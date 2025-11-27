a = [32.0	2.9341536800865753 33.0	10.750000000000009 100.0	10.750000000000002 3.0	16.989109302995107 57.14364771409426	36.050000000000004 87.0	25.450000000000006 3.8186624191061673	49.6 79.0	49.10534630385387 70.65	49.1];	
b = reshape(a, 2,[]);

mass = b(1,:) / 1000;  % g to kg
position = b(2,:) / 100; % cm to meters 

iyyzz = sum(mass.*position.^2); % about the nose of the rocket 

Icg = iyyzz - sum(mass) * .28^2; 
% ^^ Parallel axis theorem ^^ % 
