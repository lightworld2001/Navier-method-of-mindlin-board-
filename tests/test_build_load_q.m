% Test script for build_load_q
% Run with: octave -qf tests/test_build_load_q.m

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repo_root);

% Shared inputs
a = 1;
b = 1;
m_vec = [1; 2];
n_vec = [1; -1];
num = numel(m_vec);
syms x y real

% Scenario 1: constant
opts = struct('mode', 'constant', 'q0', 3.5);
load_const = build_load_q(a, b, m_vec, n_vec, num, x, y, opts);
assert(isequal(load_const.q_expr, opts.q0), 'constant: q_expr mismatch');
assert(isempty(load_const.Q_vals), 'constant: Q_vals should be empty');

% Scenario 2: fourier_numeric
Q_vals = [2; -1];
opts = struct('mode', 'fourier_numeric', 'Q_vals', Q_vals);
load_fourier = build_load_q(a, b, m_vec, n_vec, num, x, y, opts);
assert(numel(load_fourier.Q_vals) == num, 'fourier_numeric: Q_vals length mismatch');
assert(~isempty(load_fourier.q_expr), 'fourier_numeric: q_expr should be non-empty');

% Optional Scenario 3: projection
opts = struct('mode', 'projection', 'q0', 5, 'region', [0, a, 0, b]);
load_proj = build_load_q(a, b, m_vec, n_vec, num, x, y, opts);
assert(~any(isnan(load_proj.Q_vals)), 'projection: Q_vals contains NaN');

fprintf('All build_load_q tests passed.\n');
