% Parameters
numNuclei = 10;             % Number of nuclei
Zl = 8; Zr = 12;            % Coordinates of the end of NMJ
cellLength = 20;            % Length of the cell
numSteps = 1000;            % Number of time steps
deltaT = 0.1;               % Time step size
repulsionMinStrength = 1; % Minimum strength of repulsion
repulsionMaxStrength = 10;   % Maximum strength of repulsion
diffusionCoefficient = 2;   % Diffusion coefficient for signal molecules
nucleiDiffusionCoeff = 0.02; % Diffusion coefficient for nuclei
L = 1;                      % Range of force
Secretion = 50;             % Secretion rate of the signal molecules
Degrade = 0.01;             % Degradation rate of signal molecules
initialNumSignalMolecules = 20;

% loop for simulating a few cells
for c=1:10

% Initialize positions randomly within the cell
nucleiPositions = cellLength * rand(numNuclei, 1)';

% Initialize vector of repulsion strength
repulsionStrength = nucleiPositions;

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
        repulsionStrength(i) = 1./(repulsionMinStrength + (repulsionMaxStrength - repulsionMinStrength) * (nearbyMolecules / length(signalPositions)));
    end

    for i = 1:numNuclei
        % Repulsion: adjust positions to maintain minimum distance
        F = sum(repulsionStrength(:) .* sign(newNucleiPositions(i) - newNucleiPositions(:)) .* exp(-abs(newNucleiPositions(i) - newNucleiPositions(:)) / L));
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
%figure;
%hold on;
%for i = 1:numNuclei
%    plot(1:numSteps, nucleiPositionHistory(:, i));
%end
%xlabel('Time Steps');
%ylabel('Position');
%title('Random Diffusion of Nuclei in 1D Cell with Signal-Dependent Repulsion');
%legend(arrayfun(@(x) ['Nucleus ' num2str(x)], 1:numNuclei, 'UniformOutput', false));
%hold off;

% Recording ALL nuclear positions over the whole time
pos = [];
for k = 500:numSteps
    pos = [pos, nucleiPositionHistory(k, :)];
end

dx=cellLength/20;x=(1:dx:cellLength); density=x;
for j=1:length(x) density(j)=length(find(abs(pos-x(j))<dx)); end

%subplot(2,1,1)
%hist(pos) 
%subplot(2,1,2)
%hist(signalPositions)

%plot(newNucleiPositions,repulsionStrength(:),'.')
%density=density/(mean(density)); den=(1./density)/(mean((1./density)));
%plot(x(2:19),density(2:19),x(2:19),den(2:19))

dist=(1:numNuclei);
x=sort(newNucleiPositions); for j=2:numNuclei-1 dist(j) = (x(j+1)-x(j-1))/2; end
dist(1) = (x(1)+x(2))/2; 
dist(numNuclei) = cellLength - (x(numNuclei)+x(numNuclei-1))/2; dist=dist/mean(dist);
xx = x/19 - (0.5+1/19); %plot(xx,dist,'*'), axis([-0.5 0.5 0 2])

xxx=[xxx xx];ddist=[ddist dist]; % adding results for many cells

end

mm = fitlm(xxx,ddist,[0;1;2;3]);
estCoeffs4=mm.Coefficients.Estimate;pVals4=mm.Coefficients.pValue;
plot(xxx,ddist,'*k',xxx, estCoeffs4(1)+estCoeffs4(2)*xxx+estCoeffs4(4)*xxx.^3,'ob'), axis([-0.5 0.5 0 2])
pp4=[pVals4(1),pVals4(2),pVals4(3),pVals4(4)], estCoeffs4(2)