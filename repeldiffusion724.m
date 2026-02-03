% Parameters
numNuclei = 20;             % Number of nuclei
cellLength = 20;            % Length of the cell
numSteps = 800;            % Number of time steps
deltaT = 0.1;               % Time step size
repulsionMinStrength = 0.5; % Minimum strength of repulsion
repulsionMaxStrength = 2;   % Maximum strength of repulsion
diffusionCoefficient = 0.2; % Diffusion coefficient
L = 1;                      % Range of force

% Initialize positions randomly within the cell
positions = cellLength * rand(numNuclei, 1)';

% Store positions for plotting
positionHistory = zeros(numSteps, numNuclei);
positionHistory(1, :) = positions;

% Simulation loop
for step = 2:numSteps
    % Calculate random displacements
    displacements = sqrt(2 * diffusionCoefficient * deltaT) * randn(numNuclei, 1);
    
    % Update positions with displacements
    newPositions = positions + displacements';
    
    % Enforce boundary conditions (nuclei cannot leave the cell)
    newPositions = max(newPositions, 1);
    newPositions = min(newPositions, cellLength);
    
    % Calculate the current repulsion strength as a sine function of time
    repulsionStrength = repulsionMinStrength + (repulsionMaxStrength - repulsionMinStrength) *... 
        (sin(pi * (-1+2*newPositions(:)/cellLength)) + 1) / 2;
    
    % Repulsion: adjust positions to maintain minimum distance
    for i = 1:numNuclei
        F = sum(repulsionStrength .* sign(newPositions(i) - newPositions(:)) .* exp(-abs(newPositions(i) - newPositions(:)) / L));
        newPositions(i) = newPositions(i) + deltaT * F;
    end
    
    % Enforce boundary conditions again after repulsion adjustments
    newPositions = max(newPositions, 1);
    newPositions = min(newPositions, cellLength);
    
    % Update positions
    positions = newPositions;
    
    % Store positions for plotting
    positionHistory(step, :) = positions;
end

% Plot positions
figure;
hold on;
for i = 1:numNuclei
    plot(1:numSteps, positionHistory(:, i));
end
xlabel('Time Steps');
ylabel('Position');
title('Random Diffusion of Nuclei in 1D Cell with Sine Wave Repulsion');
legend(arrayfun(@(x) ['Nucleus ' num2str(x)], 1:numNuclei, 'UniformOutput', false));
hold off;

% Recording ALL nuclear positions over the whole time
pos = [];
for k = 500:numSteps
    pos = [pos, positionHistory(k, :)];
end