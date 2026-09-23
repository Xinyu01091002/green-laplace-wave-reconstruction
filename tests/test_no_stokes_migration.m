function test_no_stokes_migration()
% Guard the public runtime against reintroducing diagonal repairs.
root=setup_green_laplace();
files=dir(fullfile(root,'src','internal','*.m'));
for k=1:numel(files)
    body=fileread(fullfile(files(k).folder,files(k).name));
    forbidden={'pair_stokes_corrections(', ...
        'third_order_stokes_correction(', ...
        'third_order_surface_stokes_correction(', ...
        'apply_stokes_corrections','validation/oracles','SWORD-VWA'};
    assert(~any(contains(body,forbidden)),files(k).name);
end
n=32;mode=[0:n/2-1,-n/2:-1];[qx,qy]=meshgrid(mode,mode);
s=complex(zeros(n));s(1,3)=n^2*0.01;s(2,4)=n^2*0.006*exp(0.3i);
[eta,a,parts,psi]=gl_no_stokes_eta33(s,qx,qy,2,4,root,4);
assert(~a.stokes_correction_used && strcmp(a.crossing_gate,'none'));
assert(norm(parts.crossing_correction_spectrum(:))==0);
assert(norm(parts.Psi33_crossing_correction_spectrum(:))==0);
assert(all(isfinite(eta),'all') && all(isfinite(psi),'all'));
[lower,lower_a]=finite_depth_directional_nested_green_laplace_eta33(s,qx,qy,2,root);
assert(norm(eta(:)-lower(:))/norm(lower(:))<1e-10);
assert(~lower_a.stokes_correction_used);
[e4,p4,a4]=gl_pure_sum_order4(s,qx,qy,2,4);
assert(all(isfinite(e4),'all') && all(isfinite(p4),'all'));
assert(~a4.order4_stokes_correction_used && ~a4.all_nested_stokes_corrections_used);
assert(norm(e4(:))>0 && norm(p4(:))>0);
end
