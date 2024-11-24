clc; clear;

data = readmatrix("results/alpha-0.5_dletat-1e-2\output_u.txt");
meta_data = readtable("results/alpha-0.5_dletat-1e-2\mata_data.txt");
y = linspace(meta_data.y_0, meta_data.y_N, meta_data.N);

for i = 1:length(data(:,1)) 
    plot(y, linspace(0,1,length(y)), "Color","k");
    hold on
    ylim([-1,1]);
    plot(y, data(i, :), "LineWidth",2, "Color","b");
    hold off
    title("$\Delta t$", "Interpreter","latex")

    drawnow
    input("press");
end


