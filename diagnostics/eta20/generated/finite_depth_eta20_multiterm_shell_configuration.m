function config = finite_depth_eta20_multiterm_shell_configuration()
% Generated only by the Wolfram multi-term/shell freeze.
config = struct();
config.endpoint_first_exponents = [-1 1 0 -1 3 2 1 0 -1];
config.endpoint_second_exponents = [-1 -1 0 1 -1 0 1 2 3];
config.endpoint_coefficients = [-1.5 -0.375 0.125 -0.375 -0.053472222222222222 -0.088888888888888889 0.20625 -0.088888888888888889 -0.053472222222222222];
config.shell_edges = [0 0.1 0.3 Inf];
config.shell_anchor_q = [0.05 0.2];
config.shell_orders = [8 6 4];
config.slow_balance_channels = [12 8 1];
config.fast_balance_channels = [1 1 1];
config.multiterm_h42_transform_count = 1341;
config.multiterm_h42_product_count = 2897;
config.shell864_transform_count = 12445;
config.shell864_product_count = 59213;
config.freeze_sha256 = '1fb4670e0b7966becaa10540b4798deae6dcf84e2c88071853b8f93e01dcb936';
config.overall_exact_gate_pass = true;
end
