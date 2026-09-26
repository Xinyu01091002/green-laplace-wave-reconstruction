function summarize_audit(out)
load(fullfile(out,'audit.mat'),'T','cases');angular={};temporal={};
for p=1:5
 for pair=[15 7.5;7.5 3.75;3.75 1.875].'
  a=cases{find(T.probe==p&T.N==1101&T.angle_deg==pair(1)&T.taper_s==0,1)};
  b=cases{find(T.probe==p&T.N==1101&T.angle_deg==pair(2)&T.taper_s==0,1)};
  angular(end+1,:)={p,pair(1),pair(2),norm(a.pred(b.mask)-b.pred(b.mask))/norm(b.pred(b.mask))};
 end
end
for taper=[0 10 20]
 a=cases{find(T.probe==1&T.N==1101&T.angle_deg==1.875&T.taper_s==taper&~T.full_Hilbert_prefix,1)};
 for n=[1051 1001]
  b=cases{find(T.probe==1&T.N==n&T.angle_deg==1.875&T.taper_s==taper&~T.full_Hilbert_prefix,1)};ix=find(b.mask);
  temporal(end+1,:)={taper,n,b.t(end),norm(b.pred(ix)-a.pred(ix))/norm(a.pred(ix)),norm(b.filtered(ix,1)-a.filtered(ix,1))/norm(a.filtered(ix,1))};
 end
end
A=cell2table(angular,'VariableNames',{'probe','coarse_deg','fine_deg','prediction_change_L2'});B=cell2table(temporal,'VariableNames',{'taper_s','N','end_s','prediction_change_L2','reference_change_L2'});
writetable(A,fullfile(out,'angular_changes.csv'));writetable(B,fullfile(out,'temporal_changes.csv'));disp(A);disp(B);
end
