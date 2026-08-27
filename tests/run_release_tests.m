function summary = run_release_tests()
%RUN_RELEASE_TESTS Execute the self-contained Green--Laplace release checks.
project_root = setup_green_laplace();
addpath(fullfile(project_root,'tests'));
addpath(fullfile(project_root,'diagnostics','eta20'));
addpath(fullfile(project_root,'diagnostics','eta20','generated'));
tests = {@test_linear_frontend,@test_eta22_ordered_parity, ...
    @test_unified_order3_surface,@test_domain_guard, ...
    @test_frozen_interfaces,@test_eta20_diagnostics};
records = repmat(struct('name','','pass',false,'message',''),numel(tests),1);
for index = 1:numel(tests)
    records(index).name = func2str(tests{index});
    try
        tests{index}();
        records(index).pass = true;
        records(index).message = 'pass';
        fprintf('PASS %s\n',records(index).name);
    catch exception
        records(index).message = exception.getReport('basic','hyperlinks','off');
        fprintf(2,'FAIL %s\n%s\n',records(index).name,records(index).message);
    end
end
summary = struct('matlab_release',version('-release'), ...
    'records',records,'overall_pass',all([records.pass]));
if ~summary.overall_pass
    error('green_laplace:TestsFailed','One or more release tests failed.');
end
end

function test_linear_frontend()
c = gl_spectral_coefficients(1,9.81,1,0.01,0,2,0,0,0);
[eta,psi,~,~,parts,audit] = gl_spectral_surface(c,2*pi,2*pi,32,32,0);
x = (0:31)*(2*pi/32);
expected_eta = repmat(0.01*cos(2*x),32,1);
assert(relative_error(eta,expected_eta)<5e-14);
expected_psi = repmat(real(-1i*9.81/c.intrinsic_omega*0.01 ...
    .*exp(1i*2*x)),32,1);
assert(relative_error(psi,expected_psi)<5e-14);
assert(isfield(parts,'psi11') && ~audit.flat_potential_exposed);
end

function test_eta22_ordered_parity()
n = 64;
mode = [0:(n/2-1),-n/2:-1];
[mx,my] = meshgrid(mode,mode);
qx = 0.25*mx;
qy = 0.25*my;
spectrum = complex(zeros(n));
modes = [3,0;4,1;5,-1;6,2];
amplitudes = [0.017*exp(0.13i);0.012*exp(-0.41i); ...
    0.009*exp(0.77i);0.006*exp(-1.04i)];
for index = 1:size(modes,1)
    spectrum(mod(modes(index,2),n)+1,mod(modes(index,1),n)+1) = ...
        n*n*amplitudes(index);
end
root = string(fileparts(fileparts(mfilename('fullpath'))));
[fast,~] = green_laplace_eta22(spectrum,qx,qy,1,8,root);
direct = finite_depth_directional_pure_gl8_ordered_pair( ...
    spectrum,qx,qy,1,root);
assert(relative_error(fast,direct)<3e-11);
end

function test_unified_order3_surface()
c = gl_spectral_coefficients(3,9.81,1,[0.01,0.006],[0,0.002], ...
    [2,3],[0,1],0,0,struct('eta22_rank',4));
[eta,psi,~,~,parts,audit] = gl_spectral_surface( ...
    c,2*pi,2*pi,64,64,0);
assert(all(isfinite(eta),'all') && all(isfinite(psi),'all'));
required = {'eta22','psi22','eta33','psi33'};
assert(all(isfield(parts,required)));
assert(all(ismember(required,audit.included_components)));
assert(~isfield(parts,'Phi33') && ~audit.flat_potential_exposed);
end

function test_domain_guard()
failed = false;
try
    gl_spectral_coefficients(3,9.81,1,0.01,0,0.4,0,0,0);
catch exception
    failed = strcmp(exception.identifier,'green_laplace:Order3Domain');
end
assert(failed);
end

function test_frozen_interfaces()
root = string(fileparts(fileparts(mfilename('fullpath'))));
names = [ ...
    "finite_depth_directional_order2_eta22_pure_gl8.json", ...
    "finite_depth_directional_order2_psi22_dual_branch_gl.json", ...
    "finite_depth_directional_order3_crossing_stokes_gate.json", ...
    "finite_depth_directional_order3_nested_green_laplace.json"];
for index = 1:numel(names)
    value = jsondecode(fileread(fullfile( ...
        root,'symbolic','generated',names(index))));
    assert(value.overall_exact_gate_pass);
    if isfield(value,'candidate_id')
        assert(startsWith(lower(string(value.candidate_id)),"gl-"));
    end
end
eta20=jsondecode(fileread(fullfile( ...
    root,'symbolic','generated','eta20_neumann_r_series.json')));
assert(eta20.overall_exact_gate_pass && ~eta20.oracle_or_mf12_used ...
    && ~eta20.sampled_selection_used);
for power=[2,4,6]
    record=eta20.records.(sprintf('x%d',power));
    assert(startsWith(string(record.candidate_id),"neumann-eta20-r"));
    assert(record.neumann_max_power==power ...
        && record.neumann_layers==power+1 ...
        && record.inverse_identity_pass && record.residual_identity_pass ...
        && record.compiler_identity_pass);
end
end

function test_eta20_diagnostics()
n=64;dk=0.25;mode=[0:(n/2-1),-n/2:-1];
[mx,my]=meshgrid(mode,mode);kx=dk*mx;ky=dk*my;
spectrum=complex(zeros(n));
modes=[3,0;4,1;5,-1;6,2];
amplitudes=[0.017*exp(0.13i);0.012*exp(-0.41i); ...
    0.009*exp(0.77i);0.006*exp(-1.04i)];
for index=1:size(modes,1)
    row=mod(modes(index,2),n)+1;
    column=mod(modes(index,1),n)+1;
    mirror_row=mod(-modes(index,2),n)+1;
    mirror_column=mod(-modes(index,1),n)+1;
    spectrum(row,column)=n*n*amplitudes(index)/2;
    spectrum(mirror_row,mirror_column)=conj(spectrum(row,column));
end
root=string(fileparts(fileparts(mfilename('fullpath'))));
for rank=[6,12,16]
    [field,audit]=eta20_green_laplace_shared( ...
        spectrum,kx,ky,1,"shared"+string(rank),root);
    assert(all(isfinite(field),'all'));
    assert(strcmp(audit.information_boundary,'eta11-only'));
end
for power=[2,4,6]
    [field,audit,state]=eta20_neumann_r_series_ordered_pair( ...
        spectrum,kx,ky,1,power,root);
    assert(all(isfinite(field),'all') && isreal(field));
    assert(~audit.production_candidate && ~audit.fixed_fft_validated);
    assert(state.eta20_spectrum(1,1)==0);
end
diagnostics=gl_supported_diagnostics();
assert(height(diagnostics)==6 && ~any(diagnostics.FixedFFTValidated(4:6)));
end

function value = relative_error(candidate,reference)
value = norm(candidate(:)-reference(:))/max(norm(reference(:)),realmin);
end
