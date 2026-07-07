clear; clf; hold on;
A = load('./Heat1D_MG.m');
Nx = 128; h = 1/Nx;
x = linspace(0.5*h, 1 - 0.5*h, Nx);
m = length(A)/Nx;
for i = 1:m
    plot(x, A(1+(i-1)*Nx:i*Nx), 'k', 'linewidth', 1);
end

set(gca, 'xtick', [0 0.5 1]); set(gca, 'ytick', [-1 0 1]);
set(gca, 'fontsize', 25); axis([0 1 -1 1]);
text('Interpreter', 'latex', 'String', '$phi$', 'Position', [-0.07 0.7], 'FontSize', 30)
text('Interpreter', 'latex', 'String', '$x$', 'Position', [0.75 -1.15], 'FontSize', 30)