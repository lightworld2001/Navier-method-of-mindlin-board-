clc
clear

% Mindlin板能量变分法-Navier级数优化解（含转角）
% 简支矩形板，均布荷载

% 板参数
a = 1;      % 板长(m)
b = 1;      % 板宽(m)
h = 0.02;   % 板厚(m)
E = 2.1e11; % 弹性模量(Pa)
nu = 0.3;   % 泊松比
q = 50000;   % 均布荷载(N/m^2)
kappa = 5/6;% 剪切修正系数

D = E*h^3/(12*(1-nu^2));
G = E/(2*(1+nu));

% 级数项数与下标范围
M = 10; N = 10;
m_min = 1; n_min = 1;
m_max = M; n_max = N;
mn_list = [];
for m = m_min:1:m_max
    if m == 0
        continue;
    end
    for n = n_min:1:n_max
        if n == 0
            continue;
        end
        mn_list = [mn_list; m n];
    end
end
num = size(mn_list,1);

% 优化变量：A为num行3列矩阵，每列分别为A_w, A_tx, A_ty
A0 = zeros(num,3); % 初始值

% 优化求解（增加nu参数传递）
energy_fun = @(A) energyfunctional(A, mn_list, a, b, h, D, G, kappa, q, nu);
grad_fun = @(A) energyfunctional_grad(reshape(A,num,3), mn_list, a, b, h, D, G, kappa, q, nu);
options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10,'MaxIter',1000);
A = fsolve(grad_fun, A0(:), options); % A为列向量，reshape为矩阵
A = reshape(A,num,3);

% 拆分系数
A_w = A(:,1);
A_tx = A(:,2);
A_ty = A(:,3);

% 计算场分布
Nx = 50; Ny = 50;
x = linspace(0,a,Nx);
y = linspace(0,b,Ny);
[X,Y] = meshgrid(x,y);
w = zeros(Ny,Nx); tx = zeros(Ny,Nx); ty = zeros(Ny,Nx);

for k = 1:num
    m = mn_list(k,1); n = mn_list(k,2);
    if m > 0 && n > 0
        phi = sin(m*pi*X/a) .* cos(n*pi*Y/b);
    elseif m < 0 && n < 0
        phi = cos(abs(m)*pi*X/a) .* sin(abs(n)*pi*Y/b);
    elseif m > 0 && n < 0
        phi = sin(m*pi*X/a) .* sin(abs(n)*pi*Y/b);
    elseif m < 0 && n > 0
        phi = cos(abs(m)*pi*X/a) .* cos(n*pi*Y/b);
    else
        phi = zeros(size(X));
    end
    w = w + A_w(k)*phi;
    tx = tx + A_tx(k)*phi;
    ty = ty + A_ty(k)*phi;
end

% 以A0为零矩阵，梯度检验（积分形式）
A0 = zeros(num,3);
eps = 1e-6;
grad = zeros(size(A0));
for i = 1:size(A0,1)
    for j = 1:size(A0,2)
        A1 = A0;
        A1(i,j) = A1(i,j) + eps;
        grad(i,j) = (energyfunctional(A1, mn_list, a, b, h, D, G, kappa, q, nu) - ...
                     energyfunctional(A0, mn_list, a, b, h, D, G, kappa, q, nu)) / eps;
    end
end
fprintf('A=0时能量泛函梯度最大值: %.6e\n', max(abs(grad(:))));

fprintf('最大挠度 w_max = %.6e\n', max(abs(w(:))));
fprintf('最大转角 theta_x_max = %.6e\n', max(abs(tx(:))));
fprintf('最大转角 theta_y_max = %.6e\n', max(abs(ty(:))));
Pi = energyfunctional(A, mn_list, a, b, h, D, G, kappa, q, nu);
fprintf('最终总能量 Pi = %.6e\n', Pi);

% 绘图
figure;
subplot(1,3,1); surf(X,Y,w); title('挠度w'); xlabel('x'); ylabel('y'); zlabel('w'); colorbar;
subplot(1,3,2); surf(X,Y,tx); title('\theta_x'); xlabel('x'); ylabel('y'); zlabel('\theta_x'); colorbar;
subplot(1,3,3); surf(X,Y,ty); title('\theta_y'); xlabel('x'); ylabel('y'); zlabel('\theta_y'); colorbar;