function [a,info]=prepare_directional_joint_wavegroup( ...
        initial,first,t,s,probe,width,target_energy,minimum_count)
%PREPARE_DIRECTIONAL_JOINT_WAVEGROUP Peak-centred observed wave-group band.
arguments
    initial (1,1) struct
    first (:,1) double
    t (:,1) double
    s (1,1) struct
    probe (1,1) double {mustBeInteger,mustBePositive}
    width (1,1) double {mustBePositive} = 7.5
    target_energy (1,1) double {mustBeGreaterThan(target_energy,0),mustBeLessThanOrEqual(target_energy,1)} = .999
    minimum_count (1,1) double {mustBeInteger,mustBePositive} = 20
end

N=numel(t);dt=mean(diff(t));tr=t-t(1);F=fft(first)/N;
use=abs(initial.C)>0 & initial.kx>0;
C=initial.C(use);kx0=initial.kx(use);ky0=initial.ky(use);w0=initial.om(use);
C=C(:);kx0=kx0(:);ky0=ky0(:);w0=w0(:);
all_bins=(1:floor((N-1)/2))';
positive_energy=2*abs(F(all_bins+1)).^2;
total_energy=sum(abs(F).^2);
[~,peak_index]=max(positive_energy);
lower=peak_index;upper=peak_index;
while sum(positive_energy(lower:upper))/total_energy<target_energy ...
        || upper-lower+1<minimum_count
    left_energy=-Inf;right_energy=-Inf;
    if lower>1,left_energy=positive_energy(lower-1);end
    if upper<numel(all_bins),right_energy=positive_energy(upper+1);end
    if ~isfinite(left_energy) && ~isfinite(right_energy)
        error('The requested wave-group energy/count cannot be reached.');
    elseif left_energy>=right_energy
        lower=lower-1;
    else
        upper=upper+1;
    end
end
bins=all_bins(lower:upper);
w=2*pi*bins/(N*dt);
observed=2*conj(F(bins+1));
input=real(exp(-1i*tr*w.')*observed);
k=arrayfun(@(v)fzero(@(q)s.g*q*tanh(q*s.h)-v^2, ...
    [0,max(1,2*v^2/s.g+2/s.h)]),w);

xy=s.probes_xy_m(probe,:);
A0=C.*exp(1i*(kx0*xy(1)+ky0*xy(2)));
angles=atan2(ky0,kx0)*180/pi;
centers=(-90+width/2:width:90-width/2)';
group=1+floor((angles+90)/width);
prior_time=complex(zeros(N,numel(centers)));
for direction=1:numel(centers)
    active=find(group==direction);
    for start=1:256:numel(active)
        ids=active(start:min(start+255,numel(active)));
        prior_time(:,direction)=prior_time(:,direction) ...
            +exp(-1i*tr*w0(ids).')*A0(ids);
    end
end
prior=ifft(prior_time,[],1);
[joint,condition]=allocate_directional_record(prior(bins+1,:),observed);
[K,TH]=ndgrid(k,deg2rad(centers));
B=repmat(bins,1,numel(centers));
KX=K.*cos(TH);KY=K.*sin(TH);
active=joint~=0;
a=struct('A',joint(active),'omega',2*pi*B(active)/(N*dt), ...
    'kx',KX(active),'ky',KY(active));

retained_energy=sum(positive_energy(bins))/total_energy;
projection_relative=norm(input-first)/norm(first);
info=struct( ...
    'condition',condition,'input',input,'retained_bins',bins, ...
    'frequency_count',numel(bins),'peak_bin',all_bins(peak_index), ...
    'band_rad_s',[min(w),max(w)], ...
    'initial_band_rad_s',[min(w0),max(w0)], ...
    'projection_relative',projection_relative,'angle_deg',width, ...
    'target_energy_fraction',target_energy, ...
    'retained_energy_fraction',retained_energy, ...
    'excluded_DC_energy_fraction',abs(F(1))^2/total_energy, ...
    'selection_rule','peak-centred contiguous band; greedily add stronger adjacent bin; enforce energy and minimum count', ...
    'minimum_frequency_count',minimum_count, ...
    'min_parent_kh',min(k*s.h),'max_parent_kh',max(k*s.h));
if retained_energy<target_energy || numel(bins)<minimum_count
    error('Wave-group input selection did not satisfy its declared gates.');
end
if abs(projection_relative^2-(1-retained_energy))>1e-10
    error('Observed first-harmonic Parseval audit failed.');
end
end
