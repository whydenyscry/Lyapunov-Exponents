function [t, xsol, lyap_exp] = odeExplicitGeneralLE(c_vector, A_matrix, b_vector, ode_fun, jacobian_fun , tspan, tau, incond)
s_stages = length(c_vector);
m = length(incond);

c_vector = reshape(c_vector, [s_stages 1]);
b_vector = reshape(b_vector, [s_stages 1]);
incond = reshape(incond, [m 1]);

t = (tspan(1) : tau : tspan(2)) .';
xsol = zeros(m + m ^ 2, length(t));
xsol(:, 1) = [incond; reshape(eye(m), [m ^ 2 1])];

lyap_sum = zeros(m, 1);
lyap_exp = zeros(m, length(t));

for n = 1:length(t)-1
    t_n = t(n);
    x_n = xsol(:, n);

    K_matrix = zeros(m + m ^ 2, s_stages);
    K_matrix(1 : m, 1) = ode_fun(t_n, x_n(1 : m));
    K_matrix(m + 1 : end, 1) = reshape(jacobian_fun(t_n, x_n(1 : m)) * reshape(x_n(m + 1 : end), [m m]), [m ^ 2 1]);

    for i = 2:s_stages
        t_arg = t_n + tau * c_vector(i);
        x_arg = x_n + tau * K_matrix(:, 1 : i - 1) * A_matrix(i, 1 : i - 1) .';
        K_matrix(1 : m, i) = ode_fun(t_arg, x_arg);
        K_matrix(m + 1 : end, i) = reshape(jacobian_fun(t_arg, x_arg(1 : m)) * reshape(x_arg(m + 1 : end), [m m]), [m ^ 2 1]);
    end

    xsol(:, n + 1) = x_n + tau * K_matrix * b_vector;

    [Q, R] = qr(reshape(xsol(m + 1 : end, n + 1), [m m]));

    xsol(m + 1 : end, n + 1) = reshape(Q, [m ^ 2, 1]);

    lyap_sum = lyap_sum + log(abs(diag(R)));

    lyap_exp(:, n + 1) = lyap_sum / (t(n + 1) - t(1));
end

LEs = sort(lyap_exp(:, end), 'descend');

if sum(LEs > 0) == m
    disp("All LEs are positive");
elseif sum(LEs < 0) == m
    disp("All LEs are negative");
else
    cum_sum_LEs = cumsum(LEs);
    for k = 1 : m - 1
        if cum_sum_LEs(k) * cum_sum_LEs(k + 1) < 0
            fprintf("D_KY ≈ %.5f\n", k + cum_sum_LEs(k) / abs(LEs(k + 1)));
        end
    end
end

xsol = xsol .';
lyap_exp = lyap_exp .';
end