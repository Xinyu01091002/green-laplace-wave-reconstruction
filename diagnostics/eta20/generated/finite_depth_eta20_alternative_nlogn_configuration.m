function config = finite_depth_eta20_alternative_nlogn_configuration()
% Generated only by the Wolfram alternative-NlogN freeze.
config = struct();
config.parametrix_low_q = 0.3;
config.parametrix_high_q = 0.6;
config.log_sinc_nodes = [0.0625 0.125 0.25 0.5 1. 2. 4. 8. 16. 32.];
config.log_sinc_weights = [0.043321698784996582 0.086643397569993164 0.17328679513998633 0.34657359027997265 0.69314718055994531 1.3862943611198906 2.7725887222397812 5.5451774444795625 11.090354888959125 22.18070977791825];
config.log_sinc_radial_balance_channels = 8;
config.log_sinc_transform_count = 8395;
config.log_sinc_product_count = 17921;
config.freeze_sha256 = '832b4368c690b5deae3d39767775c230cf2e4f24aa6274bb29995dfe5f009f8d';
config.overall_exact_gate_pass = true;
config.mf12_or_oracle_used = false;
end
