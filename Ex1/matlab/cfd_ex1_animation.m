clc; clear;

methods = ["Roe_first", "Roe_second", "MacCormack", "Beam_and_Warming"];
limiters = ["no_limiter", "van_Albada", "superbee", "van_Leer", "minmod"];

method = methods(2);
limiter = limiters(2);
u0 = 1;
u1 = 0.5;
CFL = 0.8;
b = 1;
c = 0.0;
mu = 0.0;
delta_time = 0;
theta = 0;
w = 0;

[x, data, mata_data, data_of_iter] = read_data(method, u0, u1, mu, limiter, CFL, delta_time, theta, w);

ellapsed_time = 0;

fig1 = figure ('Position',[100 150 1200 500]);
for i = 1:2:length(data(2:end,1)) 
    clf(fig1);
    ellapsed_time = ellapsed_time + data_of_iter.delta_time(i);
    subplot(1,2,1)
    if mata_data.mu == 0
        plot(x(2:end-1),linspace(mata_data.u0, u1, length(x(2:end-1))), '--k')
        title(sprintf("t = %2.5f, iter = %d , CFL = %g , method = %s , limiter = %s", ellapsed_time, i, mata_data.CFL, method, limiter),Interpreter="none")
    else
        % plot(x(2:end-1), -c/b*(1+tanh(c*(x(2:end-1)-(mata_data.x_max - mata_data.x_min)/2)./(2*mata_data.mu))),"--", "Color", "k");
        if mata_data.method == "Beam_and_Warming"
            title(sprintf("t = %2.5f, iter = %d , mu = %g , theta = %g , w = %g , method = %s", ellapsed_time, i, mu, theta, w, method),Interpreter="none")
        else
            title(sprintf("t = %2.5f, iter = %d , mu = %g , method = %s", ellapsed_time, i, mu, method),Interpreter="none")
        end
    end
    subtitle(sprintf("norm: %0.10f", data_of_iter.norm(i)))
    grid on
    grid minor
    hold on
    ylim([0.4, 1.2])
    plot(x, [(data(i,1)+data(i,2))/2,data(i,2:end-1),(data(i,end-1)+data(i,end))/2], "LineWidth",2, "Color","b");
    % plot(x, data(i,:), "LineWidth",2, "Color","b");
    hold off
    

    drawnow
    % input("press");
    if (ellapsed_time >= 0.3)
        break
    end
    % pause((data_of_iter.delta_time(i))*1/(100000));
end
subplot(1,2,2)
semilogy(data_of_iter.No, data_of_iter.norm)
grid on
grid minor
% exportgraphics(fig1, 'B&W first.png','Resolution',600);

