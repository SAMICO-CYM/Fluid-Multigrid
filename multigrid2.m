clear;
hold on; box on;
A = load('./Heat1D_MG.m');
Nx = 64; h = 1/Nx; dt = 0.5 * h^2;
x = linspace(0.5*h, 1 - 0.5*h, Nx);
m = length(A)/Nx;
for i = 1:m
    plot(x, A(1+(i-1)*Nx:i*Nx), 'ko'); % numerical solution
    plot(x, exp(-(2*pi)^2*(i-1)*40*dt)*cos(2*pi*x), 'k-'); % exact solution
end

lgd = legend('numerical solution', 'exact solution'); % 범례
set(lgd, 'fontsize', 14, 'location', 'North');
set(gca, 'fontsize', 16); axis([0 1 -1 1]) % gca: get current axes: 현재 축 가져오기
text('Interpreter', 'latex', 'String', '$u$', 'Position', [-0.05 0.84], 'FontSize', 21)
text('Interpreter', 'latex', 'String', '$x$', 'Position', [0.92 -1.105], 'FontSize', 21)
print('-deps', 'heat1D_MG.eps');