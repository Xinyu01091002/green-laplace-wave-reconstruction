function [field,audit]=mf12_eta22_pure_sum_reference( ...
        state,mf12_source)
%MF12_ETA22_PURE_SUM_REFERENCE Spectral MF12 pure-sum eta22 field.
arguments
 state (1,1) struct
 mf12_source (1,1) string
end
source_file=fullfile(mf12_source,'mf12_spectral_coefficients.m');
if ~isfile(source_file),error('Official MF12 source not found: %s',source_file);end
addpath(mf12_source,'-begin');
resolved=string(which('mf12_spectral_coefficients'));
if ~strcmpi(char(resolved),char(source_file))
 error('Unexpected MF12 source resolved: %s',resolved);
end
n=size(state.eta11_spectrum,1);grid_count=n^2;
positive=state.kx>0;
threshold=1e-13*max(1,max(abs(state.eta11_spectrum),[],'all'));
parent_indices=find(positive&abs(state.eta11_spectrum)>threshold);
% The shared generator stores eta11/h.  MF12 receives physical eta11.
coefficients_eta=state.depth*state.eta11_spectrum(parent_indices)/grid_count;
a=2*real(coefficients_eta);b=2*imag(coefficients_eta);
parent_kx=state.kx(parent_indices);parent_ky=state.ky(parent_indices);
started=tic;opts=struct('enable_subharmonic',false);
coefficients=mf12_spectral_coefficients(2,1.0,state.depth,a,b, ...
 parent_kx,parent_ky,0.0,0.0,0,opts);
delta_kx=state.kx(1,2)-state.kx(1,1);
delta_ky=state.ky(2,1)-state.ky(1,1);
reference_spectrum=complex(zeros(n));
for parent=1:numel(parent_indices)
 column=mod(round(2*parent_kx(parent)/delta_kx),n)+1;
 row=mod(round(2*parent_ky(parent)/delta_ky),n)+1;
 value=complex(coefficients.A_2(parent),coefficients.B_2(parent)) ...
  *coefficients.G_2(parent);
 reference_spectrum(row,column)=reference_spectrum(row,column)+grid_count*value;
end
pair=0;
for first=1:numel(parent_indices)
 for second=(first+1):numel(parent_indices)
  pair=pair+1;plus_index=2*pair-1;
  column=mod(round((parent_kx(first)+parent_kx(second))/delta_kx),n)+1;
  row=mod(round((parent_ky(first)+parent_ky(second))/delta_ky),n)+1;
  value=complex(coefficients.A_npm(plus_index), ...
   coefficients.B_npm(plus_index))*coefficients.G_npm(plus_index);
  reference_spectrum(row,column)=reference_spectrum(row,column)+grid_count*value;
 end
end
q=state.depth*hypot(state.kx,state.ky);
minimum_q=state.config.minimum_output_midpoint_q;
certified=q>2*minimum_q;
low_output_energy=sum(abs(reference_spectrum(~certified)).^2,'all') ...
 /max(sum(abs(reference_spectrum).^2,'all'),realmin);
reference_spectrum(~certified)=0;
field=ifft2(reference_spectrum);
audit=struct('total_seconds',toc(started), ...
 'active_positive_parent_count',numel(parent_indices), ...
 'minimum_output_midpoint_q',minimum_q, ...
 'low_output_energy_fraction_before_gate',low_output_energy, ...
 'sector','positive pure-sum eta22','parallel_pool_used',false, ...
 'external_source',resolved,'rescaling_alignment_or_gain',false);
end
