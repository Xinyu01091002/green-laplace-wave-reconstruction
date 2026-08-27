function [eta20,audit,state] = eta20_neumann_r_series_ordered_pair( ...
        eta11_spectrum,kx,ky,depth,series_power,project_root)
%ETA20_NEUMANN_R_SERIES_ORDERED_PAIR R2/R4/R6 accuracy diagnostic.
%   This independent ordered-pair evaluator is validation code, not the
%   fixed-FFT production graph. series_power must be 2, 4, or 6.

arguments
    eta11_spectrum (:,:) {mustBeNumeric}
    kx (:,:) double
    ky (:,:) double
    depth (1,1) double {mustBePositive,mustBeFinite}
    series_power (1,1) double {mustBeInteger}
    project_root (1,1) string = string( ...
        fileparts(fileparts(fileparts(mfilename('fullpath')))))
end
if ~ismember(series_power,[2,4,6])
    error('series_power must be 2, 4, or 6.');
end

prepared=eta20_prepare_real_input(eta11_spectrum,kx,ky,depth);
freeze_path=fullfile(project_root,'symbolic','generated', ...
    'eta20_neumann_r_series.json');
freeze=jsondecode(fileread(freeze_path));
field_name=sprintf('x%d',series_power);
if ~freeze.overall_exact_gate_pass || freeze.oracle_or_mf12_used ...
        || freeze.sampled_selection_used || ~isfield(freeze.records,field_name)
    error('The frozen R-series information boundary failed.');
end
record=freeze.records.(field_name);
layers=series_power+1;
if record.neumann_max_power~=series_power ...
        || record.neumann_layers~=layers ...
        || ~record.inverse_identity_pass || ~record.residual_identity_pass ...
        || ~record.angular_endpoint_pass || ~record.radial_endpoint_pass ...
        || ~record.compiler_identity_pass
    error('The requested R-series record failed an exact gate.');
end

q0=prepared.q0;
c0=(8*q0*cosh(q0)^2*(4*q0+sinh(2*q0))) ...
    /(-1+8*q0^2+cosh(4*q0)-4*q0*sinh(4*q0));
group_factor=(q0*sech(q0)^2+tanh(q0))/(2*sqrt(q0*tanh(q0)));
if any(~isfinite([c0,group_factor])) || ~isreal([c0,group_factor]) ...
        || group_factor<=0
    error('Frozen R-series carrier functions are invalid.');
end

parents=find(prepared.u_spectrum~=0);
eta20_spectrum=complex(zeros(size(eta11_spectrum)));
evaluated_pairs=0;
for first_offset=2:numel(parents)
    first=parents(first_offset);
    for second_offset=1:first_offset-1
        second=parents(second_offset);
        if prepared.q(first)>=prepared.q(second)
            high=first;low=second;
        else
            high=second;low=first;
        end
        q1=prepared.q(high);q2=prepared.q(low);
        k1=[kx(high),ky(high)];k2=[kx(low),ky(low)];
        radial1=hypot(k1(1),k1(2));radial2=hypot(k2(1),k2(2));
        cosine12=max(-1,min(1,dot(k1,k2)/(radial1*radial2)));
        [row1,column1]=ind2sub(size(eta11_spectrum),high);
        [row2,column2]=ind2sub(size(eta11_spectrum),low);
        output_row=mod((row1-1)-(row2-1),prepared.ny)+1;
        output_column=mod((column1-1)-(column2-1),prepared.nx)+1;
        q_output=prepared.q_output(output_row,output_column);
        if q_output==0,continue,end
        nu1=sqrt(q1*tanh(q1));nu2=sqrt(q2*tanh(q2));
        sigma=nu1-nu2;
        forcing_k=1i*(q1*coth(q1)*nu1 ...
            -cosine12*q2*coth(q1)*nu1 ...
            +cosine12*q1*coth(q2)*nu2-q2*coth(q2)*nu2);
        forcing_d=-q1*tanh(q1)-q2*tanh(q2)+nu1*nu2 ...
            +cosine12*coth(q1)*coth(q2)*nu1*nu2;
        inverse_dno=1/(q_output*tanh(q_output));
        kernel=-0.5*forcing_d;
        for layer=1:(layers-1)
            kernel=kernel+inverse_dno^layer*( ...
                0.5i*sigma^(2*layer-1)*forcing_k ...
                -0.5*sigma^(2*layer)*forcing_d);
        end
        kernel=kernel+0.5i*inverse_dno^layers ...
            *sigma^(2*layers-1)*forcing_k;
        z=(q1-q2)^2/q_output^2;
        kernel=kernel+group_factor^(2*layers)*c0*z^layers;
        contribution=(kernel/depth)*prepared.u_spectrum(high) ...
            *conj(prepared.u_spectrum(low))/numel(eta11_spectrum);
        eta20_spectrum(output_row,output_column)= ...
            eta20_spectrum(output_row,output_column)+contribution;
        mirror_row=mod(-(output_row-1),prepared.ny)+1;
        mirror_column=mod(-(output_column-1),prepared.nx)+1;
        eta20_spectrum(mirror_row,mirror_column)= ...
            eta20_spectrum(mirror_row,mirror_column)+conj(contribution);
        evaluated_pairs=evaluated_pairs+1;
    end
end
eta20_spectrum(1,1)=0;
eta20_complex=ifft2(eta20_spectrum);
eta20=real(eta20_complex);
imaginary_leakage=norm(imag(eta20_complex(:))) ...
    /max(1,norm(eta20(:)));
if imaginary_leakage>3e-11
    error('R-series ordered-pair output is not Hermitian.');
end
audit=struct( ...
    'status','post-freeze-ordered-pair-accuracy-diagnostic', ...
    'candidate_id',char(record.candidate_id), ...
    'series_power',series_power,'neumann_layers',layers, ...
    'information_boundary','eta11-only', ...
    'strict_zero_mode','set-to-zero-separate-sector', ...
    'active_positive_parent_count',numel(parents), ...
    'evaluated_nonzero_unordered_pairs',evaluated_pairs, ...
    'pair_loops',2,'production_candidate',false, ...
    'fixed_fft_validated',false,'selection_used',false, ...
    'imaginary_output_leakage',imaginary_leakage);
state=struct('eta20_spectrum',eta20_spectrum);
end
