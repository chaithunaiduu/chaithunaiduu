
% Main script
T = 5;
N = 25;
[a0, an, bn] = compute_fourier_coefficients(@square_wave, T, N);
t_values = linspace(-T, T, 1000);
original_wave = square_wave(t_values, T);
fourier_approximation = fourier_series(t_values, T, N, a0, an, bn);
% Define square wave function
% Plotting the results
figure;
plot(t_values, original_wave, 'LineWidth', 2);
hold on;
plot(t_values, fourier_approximation, '--', 'LineWidth', 2);
xlabel('t');
ylabel('f(t)');
title(['Fourier series approximation of a square wave', ' (N=', num2str(N), ')']);
legend('Original square wave', ['Fourier Series Approximation (N=', num2str(N), ')']);
grid on;
hold off;

function y = square_wave(t, T)
    y = ones(size(t));
    y(mod(t, T) >= T/2) = -1;
end

% Compute Fourier coefficients
function [a0, an, bn] = compute_fourier_coefficients(f, T, N)
    t = linspace(0, T, 1000);
    y = f(t, T);
    a0 = (2/T) * trapz(t, y);
    an = zeros(1, N);
    bn = zeros(1, N);
    for n = 1:N
        cos_term = cos(2 * pi * n * t / T);
        sin_term = sin(2 * pi * n * t / T);
        an(n) = (2/T) * trapz(t, y .* cos_term);
        bn(n) = (2/T) * trapz(t, y .* sin_term);
    end
end

% Fourier series function
function y = fourier_series(t, T, N, a0, an, bn)
    y = (a0 / 2) * ones(size(t));
    for n = 1:N
        y = y + an(n) * cos(2 * pi * n * t / T) + bn(n) * sin(2 * pi * n * t / T);
    end
end

% Main script

