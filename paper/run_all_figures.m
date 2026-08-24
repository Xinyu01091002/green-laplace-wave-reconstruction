function outputs = run_all_figures(mode)
%RUN_ALL_FIGURES Regenerate self-contained or full paper figures.
arguments
    mode (1,1) string {mustBeMember(mode,["self-contained","full"])} = ...
        "self-contained"
end
project_root = setup_green_laplace();
addpath(fullfile(project_root,'paper'));
addpath(fullfile(project_root,'paper','internal'));
output_directory = fullfile(project_root,'results','paper');
if ~isfolder(output_directory), mkdir(output_directory); end
outputs = struct();
outputs.scalar = figure_01_scalar_convergence(output_directory);
outputs.inputs = figure_02_input_cases(output_directory);
outputs.unified_surface = figure_02_unified_surface(output_directory);
if mode=="full"
    mf12_root = string(getenv('MF12_ROOT'));
    if strlength(mf12_root)==0
        error('Set MF12_ROOT before running the full paper reproduction.');
    end
    setup_green_laplace("MF12Root",mf12_root);
    mf12_source = string(fileparts(which('mf12_spectral_coefficients')));
    input_root = fullfile(project_root,'results','paper_inputs');
    rank_root = fullfile(project_root,'results','paper_rank_data');
    generate_eta22_rank_data(mf12_source,input_root,rank_root,false);
    outputs.rank_ladder = figure_03_eta22_rank_ladder( ...
        rank_root,output_directory, ...
        "fig_eta22_field_rank_ladder_ce_no_output_cutoff",false);
    outputs.eta22_waveform = figure_04_eta22_mf12_waveform( ...
        input_root,fullfile(output_directory,'eta22_waveform'), ...
        output_directory);
    psi33_data_root = string(getenv('GL_PSI33_DATA_ROOT'));
    if strlength(psi33_data_root)>0
        outputs.psi33_waveform = figure_05_psi33_mf12_waveform( ...
            psi33_data_root,output_directory,6);
    else
        outputs.psi33_waveform = struct('status', ...
            'skipped: set GL_PSI33_DATA_ROOT to the published field archive');
    end
end
fprintf('GREEN_LAPLACE_SELF_CONTAINED_FIGURES=COMPLETE\n');
end
