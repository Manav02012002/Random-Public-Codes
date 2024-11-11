% Parameters
beta = 0.95;          % Discount factor
theta = 0.7;          % Capital elasticity
gamma = 0.2;          % Adjustment cost parameter
p = 1.2;              % Price of capital
A_bar = 1.5;          % Productivity level
eA = exp(A_bar);      % Exponential of productivity level

% Capital grid
K_min = 30;
K_max = 80;
num_K = 301;
K_grid = linspace(K_min, K_max, num_K)';  % Capital grid (column vector)

% Investment grid (since It can be negative, we consider a range)
I_min = -20;
I_max = 50;
num_I = 301;
I_grid = linspace(I_min, I_max, num_I);   % Investment grid

% Initialize value function
V = zeros(num_K, 1);        % Initial guess for value function
V_new = zeros(num_K, 1);    % To store updated value function
policy_I = zeros(num_K, 1); % To store optimal investment policy

% Tolerance level and maximum iterations
tol = 1e-6;
max_iter = 1000;
iter = 0;
diff = Inf;

% Precompute profit function for each K
profit = eA * K_grid.^theta;

% Main iteration loop
while diff > tol && iter < max_iter
    iter = iter + 1;
    for i = 1:num_K
        K_t = K_grid(i);
        % Possible next capital K_{t+1} for all I_t
        K_next = K_t + I_grid;
        
        % Ensure K_next stays within bounds
        feasible = (K_next >= K_min) & (K_next <= K_max);
        K_next_feasible = K_next(feasible);
        I_feasible = I_grid(feasible);
        
        % Calculate current period utility for feasible I_t
        adjust_cost = (gamma/2) * I_feasible.^2;
        investment_cost = p * I_feasible;
        current_utility = profit(i) - investment_cost - adjust_cost;
        
        % Interpolate V at K_{t+1}
        V_interp = interp1(K_grid, V, K_next_feasible, 'linear', 'extrap');
        
        % Total value function
        total_value = current_utility + beta * V_interp;
        
        % Find the maximum total value and corresponding investment
        [V_new(i), idx] = max(total_value);
        policy_I(i) = I_feasible(idx);
    end
    
    % Check convergence
    diff = max(abs(V_new - V));
    V = V_new;
    
    % Display iteration info every 100 iterations
    if mod(iter, 100) == 0
        fprintf('Iteration %d, diff = %.8f\n', iter, diff);
    end
end

% Plotting the value function
figure;
plot(K_grid, V, 'LineWidth', 2);
title('Value Function V(K_t)');
xlabel('Capital Stock K_t');
ylabel('Value V(K_t)');
grid on;

% Plotting the optimal investment policy function
figure;
plot(K_grid, policy_I, 'LineWidth', 2);
title('Optimal Investment Policy Function I^*(K_t)');
xlabel('Capital Stock K_t');
ylabel('Optimal Investment I_t');
grid on;

% Display iteration summary
fprintf('Converged in %d iterations with diff = %.8f\n', iter, diff);
