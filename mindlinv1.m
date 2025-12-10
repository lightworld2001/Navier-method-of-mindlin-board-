clc; clear;

%% 参数区（全部量纲统一为mm-N/mm²）
a = 50; % mm
b = 50; % mm
h = 20;   % mm
E = 2.1e5; % N/mm²
nu = 0.3;
q = 500; % N/mm²
kappa = 5/6;
D = E*h^3/(12*(1-nu^2)); % N·mm
G = E/(2*(1+nu)); % N/mm²
M = 4; N = 4;
% 支持 m,n 从 -2 开始的语义（例如 -2:-1:...），但内部用正索引访问 mn_list
m_vals = -2:M;
n_vals = -2:N;
mn_list = [];
for mm = m_vals
    for nn = n_vals
        mn_list = [mn_list; mm nn];
    end
end
num = size(mn_list,1);
m_vec = mn_list(:,1);
n_vec = mn_list(:,2);

%% 主流程区（数值/符号推导分开）


syms x y real
A_w = sym('A_w', [num 1]);
A_tx = sym('A_tx', [num 1]);
A_ty = sym('A_ty', [num 1]);

% 试函数表达式（与能量泛函关联）
w_expr = 0; tx_expr = 0; ty_expr = 0;
% 构造基函数：支持 mn_list 中的 mm, nn（可以是负数，符号决定 sin/cos）
for idx = 1:num
    mm = m_vec(idx); nn = n_vec(idx);
    if mm > 0 && nn > 0
        phi = sin(mm*pi*x/a) * cos(nn*pi*y/b);
    elseif mm < 0 && nn < 0
        phi = cos(mm*pi*x/a) * sin(nn*pi*y/b);
    elseif mm > 0 && nn < 0
        phi = sin(mm*pi*x/a) * sin(nn*pi*y/b);
    elseif mm < 0 && nn > 0
        phi = cos(mm*pi*x/a) * cos(nn*pi*y/b);
    else
        phi = 0;
    end
    w_expr = w_expr + A_w(idx)*phi;
    tx_expr = tx_expr + A_tx(idx)*phi;
    ty_expr = ty_expr + A_ty(idx)*phi;
end

w_x = diff(w_expr, x); w_y = diff(w_expr, y);
tx_x = diff(tx_expr, x); tx_y = diff(tx_expr, y);
ty_x = diff(ty_expr, x); ty_y = diff(ty_expr, y);
U_b = 0.5*D*( (-tx_x)^2 + (-ty_y)^2 + 2*nu*(-tx_x)*(-ty_y) + (1-nu)*((-tx_y-ty_x)^2)/2 );
U_s = 0.5*kappa*G*h*( (w_x-tx_expr)^2 + (w_y-ty_expr)^2 );
W_ext = q*w_expr;




% 离散罚系数参数区（位置可自定义）
k0_norm = [0.2, 0.4, 0.6, 0.8]; % x=0侧归一化位置（0~1）
N_k0 = numel(k0_norm); % x=0处离散点数
y_k0 = k0_norm * b; % x=0处离散点位置（单位mm）
k_w0 = [1e6; 1e6; 1e6; 1e6];      % x=0侧每个点的位移刚度 N/mm
k_tx0 = [5e5; 5e5; 5e5; 5e5];     % x=0侧每个点的转角刚度tx N·mm/rad
k_ty0 = [0; 0; 0; 0];     % x=0侧每个点的转角刚度ty N·mm/rad

ka_norm = [0.2, 0.4, 0.6, 0.8]; % x=a侧归一化位置（0~1）
N_ka = numel(ka_norm); % x=a处离散点数
y_ka = ka_norm * b; % x=a处离散点位置（单位mm）
k_wa  = [1e6; 1e6; 1e6; 1e6];          % x=a侧每个点的位移刚度 N/mm
k_txa = [5e5; 5e5; 5e5; 5e5];          % x=a侧每个点的转角刚度tx N·mm/rad
k_tya = [0; 0; 0; 0];          % x=a侧每个点的转角刚度ty N·mm/rad

% 离散罚函数项（x=0和x=a两侧，每个点包含w/tx/ty三项）
penalty_sum0 = 0;
w_bdy0 = subs(w_expr, x, 0); % x=0边界
tx_bdy0 = subs(tx_expr, x, 0);
ty_bdy0 = subs(ty_expr, x, 0);
for j = 1:N_k0
    yj = y_k0(j);
    penalty_sum0 = penalty_sum0 + 0.5*k_w0(j) * subs(w_bdy0^2, y, yj) ...
        + 0.5*k_tx0(j) * subs(tx_bdy0^2, y, yj) ...
        + 0.5*k_ty0(j) * subs(ty_bdy0^2, y, yj);
end

penalty_suma = 0;
w_bdya = subs(w_expr, x, a); % x=a边界
tx_bdya = subs(tx_expr, x, a);
ty_bdya = subs(ty_expr, x, a);
for j = 1:N_ka
    yj = y_ka(j);
    penalty_suma = penalty_suma + 0.5*k_wa(j) * subs(w_bdya^2, y, yj) ...
        + 0.5*k_txa(j) * subs(tx_bdya^2, y, yj) ...
        + 0.5*k_tya(j) * subs(ty_bdya^2, y, yj);
end

boundary_energy = penalty_sum0 + penalty_suma;

Pi = int(int(U_b + U_s - W_ext, x, 0, a), y, 0, b) + boundary_energy;
eqs = [];
vars = [A_w; A_tx; A_ty];
for i = 1:length(vars)
    eqs = [eqs; diff(Pi, vars(i)) == 0];
end
% 已构造方程组（为避免占用大量内存或控制台输出，已移除自动打印/预览）。

% 求解数值系数（方程对系数是线性的，使用矩阵法更稳健）
[A_mat, B_vec] = equationsToMatrix(eqs, vars);
A_num = double(A_mat);
B_num = double(B_vec);
sol_vec = A_num \ B_num; % 解线性方程组 A*vars = B
Aw_num = sol_vec(1:num);
Atx_num = sol_vec(num+1:2*num);
Aty_num = sol_vec(2*num+1:3*num);

% 数值场表达式与绘图

Nx = 40; Ny = 40;
xv = linspace(0,a,Nx); yv = linspace(0,b,Ny);
[X,Y] = meshgrid(xv,yv);
w_num = zeros(Ny,Nx); tx_num = zeros(Ny,Nx); ty_num = zeros(Ny,Nx);
for i = 1:Nx
    for j = 1:Ny
        for k = 1:num
            mm = m_vec(k); nn = n_vec(k);
            if mm > 0 && nn > 0
                phi = sin(mm*pi*xv(i)/a) * cos(nn*pi*yv(j)/b);
            elseif mm < 0 && nn < 0
                phi = cos(mm*pi*xv(i)/a) * sin(nn*pi*yv(j)/b);
            elseif mm > 0 && nn < 0
                phi = sin(mm*pi*xv(i)/a) * sin(nn*pi*yv(j)/b);
            elseif mm < 0 && nn > 0
                phi = cos(mm*pi*xv(i)/a) * cos(nn*pi*yv(j)/b);
            else
                phi = 0;
            end
            w_num(j,i) = w_num(j,i) + Aw_num(k)*phi;
            tx_num(j,i) = tx_num(j,i) + Atx_num(k)*phi;
            ty_num(j,i) = ty_num(j,i) + Aty_num(k)*phi;
        end
    end
end

figure;
subplot(1,3,1); surf(X,Y,w_num); title('挠度w'); xlabel('x'); ylabel('y'); zlabel('w'); colorbar;
subplot(1,3,2); surf(X,Y,tx_num); title('\theta_x'); xlabel('x'); ylabel('y'); zlabel('\theta_x'); colorbar;
subplot(1,3,3); surf(X,Y,ty_num); title('\theta_y'); xlabel('x'); ylabel('y'); zlabel('\theta_y'); colorbar;
