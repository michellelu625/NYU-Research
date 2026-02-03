% Parameters
numNuclei = 10;             % Number of nuclei
Zl = 6; Zr = 10;            % Coordinates of the end of NMJ
cellLength = 20;            % Length of the cell
numSteps = 1000;            % Number of time steps
deltaT = 0.1;               % Time step size
repulsionMinStrength = 0.5; % Minimum strength of repulsion
repulsionMaxStrength = 2;   % Maximum strength of repulsion
diffusionCoefficient = 5;   % Diffusion coefficient for signal molecules
nucleiDiffusionCoeff = 0.2; % Diffusion coefficient for nuclei
L = 1;                      % Range of force
Secretion = 10;             % Secretion rate of the signal molecules
Degrade = 0.01;             % Degradation rate of signal molecules
initialNumSignalMolecules = 20;

% Initialize positions randomly within the cell
nucleiPositions = cellLength * rand(numNuclei, 1)';

% Initialize random signal molecule positions
signalPositions = Zl + (Zr - Zl) * rand(initialNumSignalMolecules, 1);

% Store positions for plotting
nucleiPositionHistory = zeros(numSteps, numNuclei);
nucleiPositionHistory(1, :) = nucleiPositions;

% Simulation loop
for step = 2:numSteps
    % Calculate random displacements for nuclei
    nucleiDisplacements = sqrt(2 * nucleiDiffusionCoeff * deltaT) * randn(numNuclei, 1);
    
    % Update positions with displacements
    newNucleiPositions = nucleiPositions + nucleiDisplacements';
    
    % Enforce boundary conditions (nuclei cannot leave the cell)
    newNucleiPositions = max(newNucleiPositions, 1);
    newNucleiPositions = min(newNucleiPositions, cellLength);
    
    % Secrete new signal molecules
    newSignalPositions = Zl + (Zr - Zl) * rand(Secretion, 1);
    signalPositions = [signalPositions; newSignalPositions];
    
    % Degrade signal molecules
    aa = [];
    for j = 1:length(signalPositions)
        if rand < Degrade
            aa = [aa j];
        end
    end
    signalPositions(aa) = [];
    
    % Calculate random displacements for signal molecules
    signalDisplacements = sqrt(2 * diffusionCoefficient * deltaT) * randn(length(signalPositions), 1);
    
    % Update signal molecule positions with displacements
    signalPositions = signalPositions + signalDisplacements;
    
    % Enforce boundary conditions for signal molecules
    signalPositions = max(signalPositions, 1);
    signalPositions = min(signalPositions, cellLength);
    
    % Calculate repulsion strength based on nearby signaling molecules
    for i = 1:numNuclei
        % Find the number of nearby signaling molecules
        nearbyMolecules = sum(abs(newNucleiPositions(i) - signalPositions) < L);
        
        % Calculate the current repulsion strength
        repulsionStrength = repulsionMinStrength + (repulsionMaxStrength - repulsionMinStrength) * (nearbyMolecules / length(signalPositions));
        
        % Repulsion: adjust positions to maintain minimum distance
        F = sum(repulsionStrength .* sign(newNucleiPositions(i) - newNucleiPositions(:)) .* exp(-abs(newNucleiPositions(i) - newNucleiPositions(:)) / L));
        newNucleiPositions(i) = newNucleiPositions(i) + deltaT * F;
    end
    
    % Enforce boundary conditions again after repulsion adjustments
    newNucleiPositions = max(newNucleiPositions, 1);
    newNucleiPositions = min(newNucleiPositions, cellLength);
    
    % Update positions
    nucleiPositions = newNucleiPositions;
    
    % Store positions for plotting
    nucleiPositionHistory(step, :) = nucleiPositions;
end

% Plot nuclei positions
figure;
hold on;
for i = 1:numNuclei
    plot(1:numSteps, nucleiPositionHistory(:, i));
end
xlabel('Time Steps');
ylabel('Position');
title('Random Diffusion of Nuclei in 1D Cell with Signal-Dependent Repulsion');
legend(arrayfun(@(x) ['Nucleus ' num2str(x)], 1:numNuclei, 'UniformOutput', false));
hold off;
