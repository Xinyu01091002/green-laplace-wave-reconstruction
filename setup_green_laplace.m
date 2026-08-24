function project_root = setup_green_laplace(varargin)
%SETUP_GREEN_LAPLACE Add the public and internal MATLAB paths.
%   setup_green_laplace
%   setup_green_laplace("MF12Root",path)

parser = inputParser;
parser.addParameter("MF12Root","",@(x)ischar(x)||isstring(x));
parser.parse(varargin{:});

project_root = string(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'src'));
addpath(fullfile(project_root,'src','internal'));

mf12_root = string(parser.Results.MF12Root);
if strlength(mf12_root)>0
    candidates = [ ...
        fullfile(mf12_root,'matlab','irregularWavesMF12','Source'), ...
        fullfile(mf12_root,'irregularWavesMF12','Source'), ...
        mf12_root];
    source = "";
    for index = 1:numel(candidates)
        if isfile(fullfile(candidates(index),'mf12_spectral_coefficients.m'))
            source = candidates(index);
            break
        end
    end
    if strlength(source)==0
        error('green_laplace:MF12NotFound', ...
            'MF12Root does not contain mf12_spectral_coefficients.m.');
    end
    addpath(source,'-begin');
end
end
