clc; clear;

methods = ["Roe_first", "Roe_second", "MacCormack", "Beam_and_Warming"];
method = methods(4);
% limiter = "van_Albada";
limiter = "superbee";
% limiter = "van_Leer";
% limiter = "minmod";
% limiter = "no_limiter";
u0 = 0.0;
u1 = 1;
CFL = 0;
mu = 0.25;
delta_time = 0.5;
theta = 1;
w = 0.1;
if (mu == 0) 
    type = "inviscid";
elseif (mu ~= 0)
    type = "general";
end

if type == "inviscid"
    if method == methods(1)
        data = readmatrix(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\output_u.txt", type,CFL,u0,u1, method,CFL));
        mata_data = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\mata_data.txt", type,CFL,u0,u1, method,CFL));
        data_of_iter = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\output_iter.txt", type,CFL,u0,u1, method,CFL));
    elseif method == methods(2)
        data = readmatrix(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_%s_CFL%g\\output_u.txt", type,CFL,u0,u1, method, limiter, CFL));
        mata_data = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_%s_CFL%g\\mata_data.txt", type,CFL,u0,u1, method, limiter, CFL));
        data_of_iter = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_%s_CFL%g\\output_iter.txt", type,CFL,u0,u1, method, limiter, CFL));
    end
elseif type == "general"
    if (method == "Beam_and_Warming")
        data = readmatrix(sprintf("results\\%s\\%s\\theta%g\\w%g\\mu%g\\delta_time%g\\u0_%g_u1_%g\\output_u.txt", type, method,theta,w,mu,delta_time, u0 ,u1));
        mata_data = readtable(sprintf("results\\%s\\%s\\theta%g\\w%g\\mu%g\\delta_time%g\\u0_%g_u1_%g\\mata_data.txt", type, method,theta,w,mu,delta_time, u0 ,u1));
        data_of_iter = readtable(sprintf("results\\%s\\%s\\theta%g\\w%g\\mu%g\\delta_time%g\\u0_%g_u1_%g\\output_iter.txt", type, method,theta,w,mu,delta_time, u0 ,u1));
    else
    data = readmatrix(sprintf("results\\%s\\%s\\mu%g\\delta_time%g\\u0_%g_u1_%g\\output_u.txt", type, method,mu,delta_time, u0 ,u1));
    mata_data = readtable(sprintf("results\\%s\\%s\\mu%g\\delta_time%g\\u0_%g_u1_%g\\mata_data.txt", type, method,mu,delta_time, u0 ,u1));
    data_of_iter = readtable(sprintf("results\\%s\\%s\\mu%g\\delta_time%g\\u0_%g_u1_%g\\output_iter.txt", type, method,mu,delta_time, u0 ,u1));
    end
end
x = [mata_data.x_min, linspace(mata_data.x_min, mata_data.x_max, mata_data.N), mata_data.x_max];
ellapsed_time = 0;

fig1 = figure ('Position',[100 150 900 500]);
for i = 1:length(data(2:end,1)) 
    clf(fig1);
    ellapsed_time = ellapsed_time + data_of_iter.delta_time(i);
    subplot(1,2,1)
    if type == "inviscid"
        plot(x(2:end-1), linspace(1,mata_data.u1,length(x(2:end-1))), "Color","k");
        title(sprintf("t = %2.5f, iter = %d , CFL = %g , method = %s", ellapsed_time, i, mata_data.CFL, method),Interpreter="none")
    else
        title(sprintf("t = %2.5f, iter = %d , mu = %g , method = %s", ellapsed_time, i, mu, method),Interpreter="none")
    end
    subtitle(sprintf("norm: %0.10f", data_of_iter.norm(i)))
    grid on
    grid minor
    hold on
    plot(x, [(data(i,1)+data(i,2))/2,data(i,2:end-1),(data(i,end-1)+data(i,end))/2], "LineWidth",2, "Color","b");
    % plot(x, data(i,:), "LineWidth",2, "Color","b");
    hold off
    

    drawnow
    % input("press");
    if (ellapsed_time >= 18)
        break
    end
    pause((data_of_iter.delta_time(i))*1/(100000));
end
subplot(1,2,2)
semilogy(data_of_iter.No, data_of_iter.norm)
grid on
grid minor

