num_nuclei = 2; num_steps = 300; cell_length = 8; diffusion_coefficient = 0.1;
    %Initialize position of nuclei randomly in cell
    positions = cell_length*rand(num_nuclei, 1);

    %Time step
    dt = .1;

    %Record position at each time step 
    position_history = zeros(num_nuclei, num_steps);
    position_history(:,1) = positions;

    %Simulation loop
    for t=2:num_steps
        %Random displacement for each nucleus
        displacement = sqrt(2*diffusion_coefficient*dt)*randn(num_nuclei,1);
        %Update Positions
        positions = positions + displacement;
        %Boundary Conditions
        positions(positions<0) = -positions(positions<0);
        positions(positions>cell_length) = 2*cell_length - positions(positions>cell_length);
        position_history(:,t) = positions;
        %Distance Between
        x = sort(positions);
        distance = abs(positions(1) - positions(2));
        %Repulsion
        if distance<3
            repulsion = 2;
            x(2) = x(2) + repulsion;
            x(1) = x(1) - repulsion;
        end
        position_history(:,t) = positions;

    end

    %Plot
    figure;
    for i = 1:num_nuclei
        plot(1:num_steps, position_history(i,:));
        hold on;
    end
    xlabel('Time Step');
    ylabel('Position');
    title('Random Diffusion of Nuclei in a 1D Cell');
    hold off;