clc; clear; close all;

N = 100;
scheme = 'EXPLICIT_SW';
M_inf  = 0.9;
CFL    = 0.7;

result = read_data(N, scheme, M_inf, CFL);

fig1 = figure ('Name', '1', 'Position',[100 150 900 500]);

colors = cool(length(1:30:length(result.data)))*0.9;
j = 1;
for i = 1:N/20:length(result.data)
    i
    clf(fig1);
    hold all

    plot(result.x, result.data{i,1}.norm_u, 'Color','k')
    plot(result.x, result.data{i,1}.norm_rho, 'Color','b')
    plot(result.x, result.data{i,1}.norm_p, 'Color','r')
    j = j+1;
    
    if max(result.data{i,1}.norm_p, result.data{i,1}.norm_rho) < 3.75
        ylim([-1.3, 3.75])
    end
    xlim([result.x(1)-0.01, result.x(end)+0.01])
    grid on
    grid minor

    font_size = 15;
    title(sprintf('%s, $N = %d$, $CFL = %g$, $Re_\\infty = %g$, $M_\\infty = %g$, $t = %g$', result.metadata.scheme{1,1}, N, result.metadata.CFL, result.metadata.Re_inf, result.metadata.M_inf, result.iter_data.elapsed_time(i)),'FontSize',font_size,'Interpreter','latex')
    ylabel('ratio [-]','FontSize',font_size, "Interpreter","latex")
    xlabel('x [-]','FontSize',font_size, "Interpreter","latex")
    legend({'Velocity ratio', 'Density ratio', 'Pressure ratio'},'FontSize',font_size-4 ,'Location','southeast','Interpreter','latex')
    box on

    drawnow
    % input('press');
    % if (ellapsed_time >= 10)
    %     break
    % end
    % pause((result.iter_data.norm_delta_time(i))*1/1);
end
% exportgraphics(fig1, 'images/test.png','Resolution',400);


