function fitness = Rastrigin(x)
% Rastrigin 函数: f(x) = 10*D + Σ[x_i² - 10cos(2πx_i)],  x_i ∈ [-5.12, 5.12]
    [N, D] = size(x);
    fitness = zeros(N, 1);
    for i = 1:N
        xi = x(i, :);
        fitness(i) = 10 * D + sum(xi.^2 - 10 * cos(2 * pi * xi));
    end
end
