function [result] = read_data(N, scheme, M_inf, CFL)
    
    
    data      = readmatrix(sprintf('results\\N%d\\%s\\M%g\\CFL%g\\output_Q.txt', N, scheme, M_inf, CFL));
    metadata  = readtable(sprintf('results\\N%d\\%s\\M%g\\CFL%g\\metadata.txt', N, scheme, M_inf, CFL));
    iter_data = readtable(sprintf('results\\N%d\\%s\\M%g\\CFL%g\\output_iter_data.txt', N, scheme, M_inf, CFL));


    Qs = {};
    for i = 0:length(data(:,1))/3-1
        Qs{i+1,1} = data(i*3+1:i*3+3,2:end-1);
    end
    
    data = {};
    for i = 1:length(Qs(:,1))
        i
        data{i,1}.norm_rho = Qs{i,1}(1,:);
        data{i,1}.norm_u = Qs{i,1}(2,:) ./ Qs{i,1}(1,:);
        data{i,1}.norm_e = Qs{i,1}(3,:);
        data{i,1}.norm_p = (metadata.gamma - 1) * (Qs{i,1}(3,:) - 0.5 * Qs{i,1}(1,:) .* (Qs{i,1}(2,:) ./ Qs{i,1}(1,:)).^2);
    end

    x_min = metadata.x_min;
    x_max = metadata.x_max;
    N = metadata.N;
    L = x_max - x_min; 
    delta_x = (x_max - x_min)/(N);
    x = [];
    for i = 1:N
        x(i) = x_min + delta_x/2 + delta_x*(i-1);
    end 
    % x_new(1) = x_new(2)-delta_x;
    % x_new(end+1) = x_new(end)+delta_x;
    x = x ./ L;

    result.x = x;
    result.data = data;
    result.metadata = metadata;
    result.iter_data = iter_data;

end

