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
% 支持 m,n 从 -4 开始的语义（例如 -4:-1:...），但内部用正索引访问 mn_list
m_vals = -4:M;
n_vals = -4:N;
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

% 试函数表达式
w_expr = 0;
for idx = 1:num
    mm = m_vec(idx); nn = n_vec(idx);
    phi = basis_phi(mm, nn, x, y, a, b);
    w_expr = w_expr + A_w(idx)*phi;
end

% Kirchhoff板能量泛函
w_xx = diff(w_expr, x, 2);
w_yy = diff(w_expr, y, 2);
w_xy = diff(diff(w_expr, x), y);

U_b = 0.5*D*(w_xx^2 + w_yy^2 + 2*nu*w_xx*w_yy + 2*(1-nu)*w_xy^2);
% 将载荷 q 表示为 Fourier 级数（与基函数相同的 phi 系列），便于与试函数正交
Q = sym('Q', [num 1]); % 若有数值系数，可后续赋值 Q_vals 并做 double(Q_vals)
q_expr = 0;
for k = 1:num
    mm = m_vec(k); nn = n_vec(k);
    q_expr = q_expr + Q(k)*basis_phi(mm, nn, x, y, a, b);
end
% 指定载荷作用矩形区域（可修改），若为圆形请用极坐标数值积分
x1 = 10; x2 = 40; y1 = 10; y2 = 40;

% ---------------- 载荷输入接口（可修改） ----------------
% 符号求解：保持载荷系数 Q 为符号向量，与试函数基一一对应。
% 若需均布载荷的数值系数用于绘图，可启用如下投影计算。
compute_projection = true; % 需要时设为 true，不需要时设为 false
if compute_projection
    Q_vals = zeros(num,1);
    for k = 1:num
        mm = m_vec(k); nn = n_vec(k);
        phi_k = @(xx,yy) basis_phi(mm, nn, xx, yy, a, b);
        numInt = integral2(@(xx,yy) q .* phi_k(xx,yy), x1, x2, y1, y2);
        denInt = integral2(@(xx,yy) phi_k(xx,yy).^2, x1, x2, y1, y2);
        if denInt == 0
            Q_vals(k) = 0;
        else
            Q_vals(k) = numInt / denInt;
        end
    end
end
% ---------------- 载荷输入接口结束 ----------------

W_ext = int(int(q_expr * w_expr, x, x1, x2), y, y1, y2);

% 离散罚系数参数区（位置可自定义）
%----------------------------------------------------这里是离散罚系数参数区(刚度自定义)
k0_norm = [0.2, 0.4, 0.6, 0.8]; % x=0侧归一化位置（0~1）
N_k0 = numel(k0_norm); % x=0处离散点数
y_k0 = k0_norm * b; % x=0处离散点位置（单位mm）
k_w0 = [1e6; 1e6; 1e6; 1e6];      % x=0侧每个点的位移刚度 N/mm
k_tx0 = [5e5; 5e5; 5e5; 5e5];     % x=0侧每个点的转角刚度tx N·mm/rad

ka_norm = [0.2, 0.4, 0.6, 0.8]; % x=a侧归一化位置（0~1）
N_ka = numel(ka_norm); % x=a处离散点数
y_ka = ka_norm * b; % x=a处离散点位置（单位mm）
k_wa  = [1e6; 1e6; 1e6; 1e6];          % x=a侧每个点的位移刚度 N/mm
k_txa = [5e5; 5e5; 5e5; 5e5];          % x=a侧每个点的转角刚度tx N·mm/rad

% 离散罚函数项（x=0和x=a两侧，每个点包含w/tx/ty三项）
penalty_sum0 = 0;
w_bdy0 = subs(w_expr, x, 0); % x=0边界
w_xx_bdy0 = subs(w_xx, x, 0); % x=0边界二阶导

for j = 1:N_k0
    yj = y_k0(j);
    penalty_sum0 = penalty_sum0 ...
        + 0.5*k_w0(j) * subs(w_bdy0^2, y, yj) ...
        + 0.5*k_tx0(j) * subs(w_xx_bdy0^2, y, yj) ;
end

penalty_suma = 0;
w_bdya = subs(w_expr, x, a); % x=a边界
w_xx_bdya = subs(w_xx, x, a);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
for j = 1:N_ka
    yj = y_ka(j);
    penalty_suma = penalty_suma ...
        + 0.5*k_wa(j) * subs(w_bdya^2, y, yj) ...
        + 0.5*k_txa(j) * subs(w_xx_bdya^2, y, yj) ;
end

boundary_energy = penalty_sum0 + penalty_suma;

Pi = int(int(U_b, x, 0, a), y, 0, b) - W_ext + boundary_energy;
eqs = [];
vars = A_w;
for i = 1:length(vars)
    eqs = [eqs; diff(Pi, vars(i)) == 0];
end
[A_mat, B_vec] = equationsToMatrix(eqs, vars);
Aw_sym = A_mat \ B_vec; % 纯符号解（与 Q 等符号系数保持一致）

% 如果提供数值载荷系数 Q_vals，可用于数值绘图（不影响符号求解）
if exist('Q_vals','var') && numel(Q_vals) == num
    Aw_num = double(subs(Aw_sym, Q, Q_vals(:)));
    % 数值场表达式与绘图
    %----------------------------------------------------这里是数值场表达式与绘图
    Nx = 40; Ny = 40;
    xv = linspace(0,a,Nx); yv = linspace(0,b,Ny);
    [X,Y] = meshgrid(xv,yv);
    w_num = zeros(Ny,Nx);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:num
                mm = m_vec(k); nn = n_vec(k);
                phi = basis_phi(mm, nn, xv(i), yv(j), a, b);
                w_num(j,i) = w_num(j,i) + Aw_num(k)*phi;
            end
        end
    end

    figure;
    surf(X,Y,w_num); title('Kirchhoff板挠度w'); xlabel('x'); ylabel('y'); zlabel('w'); colorbar;
end

function phi = basis_phi(mm, nn, x, y, a, b)
    if mm > 0 && nn > 0
        phi = sin(mm*pi*x/a) .* cos(nn*pi*y/b);
    elseif mm < 0 && nn < 0
        phi = cos(mm*pi*x/a) .* sin(nn*pi*y/b);
    elseif mm > 0 && nn < 0
        phi = sin(mm*pi*x/a) .* sin(nn*pi*y/b);
    elseif mm < 0 && nn > 0
        phi = cos(mm*pi*x/a) .* cos(nn*pi*y/b);
    else
        phi = 0;
    end
end
