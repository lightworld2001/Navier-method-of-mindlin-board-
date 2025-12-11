function load = build_load_q(a,b,m_vec,n_vec,num,x,y,opts)
% BUILD_LOAD_Q  构造载荷 q(x,y) 的辅助函数（兼容 MATLAB R2016+）
% 简化实现：不使用 arguments 块以提高兼容性

if nargin < 8 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'mode'), opts.mode = 'constant'; end
if ~isfield(opts,'q0'), opts.q0 = 0; end
if ~isfield(opts,'Q_vals'), opts.Q_vals = []; end
if ~isfield(opts,'region'), opts.region = [0, a, 0, b]; end
if ~isfield(opts,'compute_projection'), opts.compute_projection = false; end

% 如果传入的 x,y 不是符号变量，重新创建符号变量 x,y
if ~isa(x,'sym') || ~isa(y,'sym')
    syms x y real
end

region = opts.region;
if numel(region) ~= 4
    error('opts.region must be [x1 x2 y1 y2]');
end
x1 = region(1); x2 = region(2); y1 = region(3); y2 = region(4);

% 内部 helper：基函数
    function phi_k = make_phi(mm, nn)
        if mm > 0 && nn > 0
            phi_k = sin(mm*pi*x/a) * cos(nn*pi*y/b);
        elseif mm < 0 && nn < 0
            phi_k = cos(mm*pi*x/a) * sin(nn*pi*y/b);
        elseif mm > 0 && nn < 0
            phi_k = sin(mm*pi*x/a) * sin(nn*pi*y/b);
        elseif mm < 0 && nn > 0
            phi_k = cos(mm*pi*x/a) * cos(nn*pi*y/b);
        else
            phi_k = sym(0);
        end
    end

% 返回结构
load = struct();
load.region = region;

mode = lower(opts.mode);
switch mode
    case 'constant'
        q0 = opts.q0;
        load.q_expr = q0;
        load.Q_vals = [];
        load.Qsym = [];
    case 'fourier_symbolic'
        Q = sym('Q', [num 1]);
        q_expr = sym(0);
        for k = 1:num
            mm = m_vec(k); nn = n_vec(k);
            q_expr = q_expr + Q(k)*make_phi(mm,nn);
        end
        load.q_expr = q_expr;
        load.Qsym = Q;
        load.Q_vals = [];
    case 'fourier_numeric'
        if isempty(opts.Q_vals) || numel(opts.Q_vals) ~= num
            error('opts.Q_vals must be a numeric vector of length num for mode=''fourier_numeric''');
        end
        Qv = double(opts.Q_vals(:));
        q_expr = sym(0);
        for k = 1:num
            mm = m_vec(k); nn = n_vec(k);
            q_expr = q_expr + Qv(k)*make_phi(mm,nn);
        end
        load.q_expr = q_expr;
        load.Q_vals = Qv;
        load.Qsym = [];
    case 'projection'
        q0 = opts.q0;
        Qv = zeros(num,1);
        for k = 1:num
            mm = m_vec(k); nn = n_vec(k);
            phi_k = make_phi(mm,nn);
            numInt = int(int(q0 * phi_k, x, x1, x2), y, y1, y2);
            denInt = int(int(phi_k^2, x, x1, x2), y, y1, y2);
            numVal = double(numInt);
            denVal = double(denInt);
            if denVal == 0 || isnan(denVal)
                Qv(k) = 0;
            else
                Qv(k) = numVal / denVal;
            end
        end
        q_expr = sym(0);
        for k = 1:num
            mm = m_vec(k); nn = n_vec(k);
            q_expr = q_expr + Qv(k)*make_phi(mm,nn);
        end
        load.q_expr = q_expr;
        load.Q_vals = Qv;
        load.Qsym = [];
    otherwise
        error('Unknown mode: %s', mode);
end

end