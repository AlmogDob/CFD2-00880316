function [x, data, mata_data, data_of_iter] = read_data(method, u0, u1, mu, limiter, CFL, delta_time, theta, w)
    if (mu == 0) 
        type = "inviscid";
    elseif (mu ~= 0)
        type = "general";
    end
    
    if type == "inviscid"
        if method == "Roe_first"
            data = readmatrix(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\output_u.txt", type,CFL,u0,u1, method,CFL));
            mata_data = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\mata_data.txt", type,CFL,u0,u1, method,CFL));
            data_of_iter = readtable(sprintf("results\\%s\\CFL%g\\u0_%g_u1_%g\\%s_CFL%g\\output_iter.txt", type,CFL,u0,u1, method,CFL));
        elseif method == "Roe_second"
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
end

