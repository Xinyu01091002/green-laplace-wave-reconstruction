function test_two_scale_migration()
root=setup_green_laplace();
addpath(fullfile(root,'research','two_scale'));
eta=run_ordered_triple_parity();
psi=run_psi33_ordered_triple_parity();
assert(eta.pass && psi.pass);
end
