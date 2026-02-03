function diffusion_simulation(num_nuclei, num_steps, cell_length, diffusion_coefficient)
    % num_nuclei: Number of nuclei
    % num_steps: Number of time steps
    % cell_length: Length of the cell
    % diffusion_coefficient: Diffusion coefficient
    
    % Initialize positions of nuclei randomly within the cell
    positions = cell_length * rand(num_nuclei, 1);
    
    % Time step (arbitrary units)
    dt = 1;
    
    % Record the positions at each time step
    positions_history = zeros(num_nuclei, num_steps);
    positions_history(:, 1) = positions;
    
    % Simulation loop
    for t = 2:num_steps
        % Calculate random displacements for each nucleus
        displacements = sqrt(2 * diffusion_coefficient * dt) * randn(num_nuclei, 1);
        
        % Update positions with periodic boundary conditions
        positions = positions + displacements;
        
        % Reflecting boundary conditions
        positions(positions < 0) = -positions(positions < 0);
        positions(positions > cell_length) = 2 * cell_length - positions(positions > cell_length);
        
        % Save positions to history
        positions_history(:, t) = positions;
    end
    
    % Plot the trajectories of the nuclei
    figure;
    for i = 1:num_nuclei
        plot(1:num_steps, positions_history(i, :));
        hold on;
    end
    xlabel('Time Step');
    ylabel('Position');
    title('Random Diffusion of Nuclei in a 1D Cell');
    hold off;
end
