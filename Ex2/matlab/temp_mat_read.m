clc; clear; close all;

data = readmatrix("results\output_Q.txt");

Qs = {};
for i = 0:length(data(:,1))/3-1
    Qs{i+1,1} = data(i*3+1:i*3+3,2:end-1);
end

data = {};
for i = 1:length(Qs(:,1))
    data{i,1}.norm_rho = Qs{i,1}(1,:);
    data{i,1}.norm_u = Qs{i,1}(2,:) ./ Qs{i,1}(1,:);
    data{i,1}.norm_e = Qs{i,1}(3,:);
end

hold all
colors = cool(length(1:30:length(data)))*0.9;
j = 1;
for i = 1:30:length(data)
    plot(data{i,1}.norm_u, "Color",colors(j,:))
    j = j+1;
end

