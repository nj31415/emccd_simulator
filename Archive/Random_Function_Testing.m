% Test script
clear;
clc;

N = 20000;
p = 0.1;

%f = @() sum(rand(N, 1) < p);
f = @() gather(sum(gpuArray.rand(N, 1) < p));
g = @() sum(rand(N, 1) < p);

trials = 100;

f_sum_time = 0;
g_sum_time = 0;

for i = 1 : trials
    f_sum_time = f_sum_time + timeit(f);
    g_sum_time = g_sum_time + timeit(g);
end

f_avg_time = f_sum_time / trials;
g_avg_time = g_sum_time / trials;

function rnd_num = binom_vec_rand(N, P)
N_form = reshape(N, 1, []);
P_form = reshape(P, size(N_form));
rnd_num = zeros(size(N_form));

parfor i = 1 : length(N_form)
    rnd_num(i) = sum(rand(N_form(i), 1) < P_form(i));
end

rnd_num = reshape(rnd_num, size(N));
end