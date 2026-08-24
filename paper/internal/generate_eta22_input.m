function state=generate_eta22_input(config)
% Configurable eta11 input for the independent eta22 playground.
% The shared eta20 routine is used only as a spectrum/discretization utility;
% no eta20 candidate, kernel, reference, or higher-order field is consumed.
arguments
 config (1,1) struct
end
base=generate_directional_input(config);
positive=base.kx>0;
analytic_spectrum=complex(zeros(size(base.eta11_spectrum)));
analytic_spectrum(positive)=2*base.eta11_spectrum(positive);
state=base;
state.eta11_analytic_spectrum=analytic_spectrum;
state.eta11_analytic=ifft2(analytic_spectrum);
state.g2=[];state.g2_audit=struct();state.g2_seconds=NaN;
state.mf12=[];state.mf12_audit=struct();
state.input_utility='shared-spectrum-discretization-only';
state.information_boundary='eta11-only';
end
