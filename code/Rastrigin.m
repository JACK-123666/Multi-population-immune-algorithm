function fitness = Rastrigin(x)
% Rastrigin 函数: f(x) = 10*D + Σ[x_i² - 10cos(2πx_i)],  x_i ∈ [-5.12, 5.12]
    D = size(x, 2);
    fitness = 10 * D + sum(x.^2 - 10 * cos(2 * pi * x), 2);
end
