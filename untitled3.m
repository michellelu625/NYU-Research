% % Parameters
numNuclei = 2;             % Number of nuclei
cellLength = 100;           % Length of the cell
numSteps = 1000;            % Number of time steps
deltaT = 0.1;               % Time step size
repulsionStrength = 1;      % Strength of repulsion
minDistance = 2;            % Minimum allowed distance between nuclei
diffusionCoefficient = 1;   % Diffusion coefficient

% Initialize positions randomly within the cell, ensuring no overlap
positions = sort(randperm(cellLength - minDistance * (numNuclei - 1), numNuclei));
positions = positions + (0:numNuclei-1) * minDistance;

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
    
    % Repulsion: adjust positions to maintain minimum distance
    for i = 1:numNuclei
        for j = 1:numNuclei
            if i ~= j
                distance = abs(newPositions(i) - newPositions(j));
                if distance < minDistance
                    % Adjust positions to maintain minimum distance
                    if newPositions(i) < newPositions(j)
                        newPositions(i) = newPositions(i) - repulsionStrength * (minDistance - distance) / 2;
                        newPositions(j) = newPositions(j) + repulsionStrength * (minDistance - distance) / 2;
                    else
                        newPositions(i) = newPositions(i) + repulsionStrength * (minDistance - distance) / 2;
                        newPositions(j) = newPositions(j) - repulsionStrength * (minDistance - distance) / 2;
                    end
                end
            end
        end
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
title('Random Diffusion of Nuclei in 1D Cell');
legend(arrayfun(@(x) ['Nucleus ' num2str(x)], 1:numNuclei, 'UniformOutput', false));
hold off;