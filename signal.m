% Parameters
Zl = 6; Zr = 10;            % Coordinates of the end of NMJ
cellLength = 20;            % Length of the cell
numSteps = 1000;            % Number of time steps
deltaT = 0.1;               % Time step size
repulsion = 1;              % Average strength of repulsion
diffusionCoefficient = 5; % Diffusion coefficient
Secretion = 10;              % Secretion of the signal rate
Degrade = 0.01;             % rate of signal molecule degradation

% Initialize random signal molecule positions
positions = Zl+(Zr - Zl) * rand(20, 1);

% Simulation loop
for step = 2:numSteps

    NewPositions = Zl+(Zr - Zl) * rand(Secretion, 1); positions = [positions; NewPositions];
    
aa=[];    for j=1:length(positions) if rand<Degrade aa=[aa j]; else aa==aa; end, end
          positions(aa)=[];

    % Calculate random displacements
    displacements = sqrt(2 * diffusionCoefficient * deltaT) * randn(length(positions),1);
    
    % Update positions with displacements
    positions = positions + displacements;
    
    % Enforce boundary conditions (nuclei cannot leave the cell)
    positions = max(positions, 1);
    positions = min(positions, cellLength);
    
 
end

% Plot positions
plot(positions,zeros(size(positions)),'o');
xlabel('x coordinate');
ylabel('Position');
title('Distribution of signaling molecules');