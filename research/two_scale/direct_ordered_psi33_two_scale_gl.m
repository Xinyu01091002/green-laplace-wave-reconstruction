function result = direct_ordered_psi33_two_scale_gl(modes,cfg)
%DIRECT_ORDERED_PSI33_TWO_SCALE_GL Independent ordered Psi33 reconstruction.

if nargin<2, cfg=struct(); end
if ~isfield(cfg,'mode'), cfg.mode='two_scale'; end
if ~isfield(cfg,'freeze_json') || isempty(cfg.freeze_json)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    cfg.freeze_json=fullfile(root,'symbolic','generated', ...
        'finite_depth_directional_order3_two_scale_gl.json');
end
freeze=jsondecode(fileread(cfg.freeze_json));
kx=modes.kx(:); ky=modes.ky(:); amplitude=modes.amplitude(:);
q=hypot(kx,ky); omega=sqrt(q.*tanh(q)); n=numel(q);

eta2=complex(zeros(n)); phi2=complex(zeros(n));
omega2=zeros(n); kx2=zeros(n); ky2=zeros(n); G2=zeros(n);
for a=1:n
    for b=1:n
        kx2(a,b)=kx(a)+kx(b); ky2(a,b)=ky(a)+ky(b);
        qout=hypot(kx2(a,b),ky2(a,b)); G2(a,b)=qout*tanh(qout);
        omega2(a,b)=omega(a)+omega(b);
        dotab=kx(a)*kx(b)+ky(a)*ky(b);
        fk=1i*(dotab/omega(a)+q(b)^2/omega(b));
        fd=-omega(b)^2+0.5*(dotab/(omega(a)*omega(b))-omega(a)*omega(b));
        r=approximate_resolvent(2,omega2(a,b),sqrt(G2(a,b)));
        eta2(a,b)=r*(G2(a,b)*fd-1i*omega2(a,b)*fk);
        phi2(a,b)=r*(-1i*omega2(a,b)*fd-fk);
    end
end

records=repmat(struct('indices',[],'kx',0,'ky',0,'coefficient',0),n^3,1);
counter=0;
for a=1:n
    for b=1:n
        for c=1:n
            counter=counter+1;
            ktotx=kx(a)+kx(b)+kx(c); ktoty=ky(a)+ky(b)+ky(c);
            qout=hypot(ktotx,ktoty); Gout=qout*tanh(qout);
            Omega=omega(a)+omega(b)+omega(c);
            Kab_dot_kc=kx2(a,b)*kx(c)+ky2(a,b)*ky(c);
            ka_dot_Kbc=kx(a)*kx2(b,c)+ky(a)*ky2(b,c);
            kb_dot_kc=kx(b)*kx(c)+ky(b)*ky(c);
            fk=1i*ka_dot_Kbc/omega(a)*eta2(b,c)-Kab_dot_kc*phi2(a,b) ...
                +1i*omega(b)*kb_dot_kc-(kx2(b,c)^2+ky2(b,c)^2)*phi2(b,c) ...
                +1i*q(c)^2/omega(c)*eta2(a,b)+0.5i*q(c)^2*omega(c);
            fd=-1i*omega2(b,c)*G2(b,c)*phi2(b,c)-omega(c)^2*eta2(a,b) ...
                -0.5*q(c)^2+1i*ka_dot_Kbc/omega(a)*phi2(b,c) ...
                -1i*omega(a)*G2(b,c)*phi2(b,c)+kb_dot_kc*omega(c)/omega(b) ...
                -omega(b)*q(c)^2/omega(c);
            r=approximate_resolvent(3,Omega,sqrt(Gout));
            phi3=r*(-1i*Omega*fd-fk);
            correction=G2(b,c)*phi2(b,c)-1i*omega(c)*eta2(a,b) ...
                -0.5i*q(c)^2/omega(c);
            kernel=phi3+correction;
            records(counter).indices=[a,b,c];
            records(counter).kx=ktotx; records(counter).ky=ktoty;
            records(counter).coefficient=kernel*amplitude(a)*amplitude(b)*amplitude(c)/4;
        end
    end
end
result=struct('records',records,'mode',cfg.mode,'target','surface_Psi33', ...
    'eta11_only',true,'external_lower_order_fields',false,'gauge','C3(t)=0');

    function r=approximate_resolvent(stage,Omega,wout)
        prefix=['order' num2str(stage)]; r=0;
        if wout==0, return; end
        use_two_scale=strcmp(cfg.mode,'two_scale') ...
            || (strcmp(cfg.mode,'inner_shared_outer_two_scale') && stage==3);
        if use_two_scale
            r=branch(freeze.primary_node_counts.(prefix).slow, ...
                number(freeze.scales.([prefix '_slow'])),1)+ ...
                branch(freeze.primary_node_counts.(prefix).fast, ...
                number(freeze.scales.([prefix '_fast'])),-1);
        else
            gl=gl_rule(freeze.shared_control_node_counts.(prefix));
            scale=number(freeze.scales.([prefix '_shared']));
            for ii=1:numel(gl.nodes)
                tau=gl.nodes(ii)/scale; aw=gl.weights(ii)*exp(gl.nodes(ii))/scale;
                r=r+aw*exp(-Omega*tau)*sinh(wout*tau)/wout;
            end
        end
        function value=branch(nn,scale,sgn)
            gl=gl_rule(nn); value=0;
            for kk=1:numel(gl.nodes)
                tau=gl.nodes(kk)/scale; aw=gl.weights(kk)*exp(gl.nodes(kk))/scale;
                value=value+sgn*aw*exp(-Omega*tau)*exp(sgn*wout*tau)/(2*wout);
            end
        end
    end
    function out=gl_rule(nn)
        source=freeze.gauss_laguerre.(['n' num2str(nn)]);
        out.nodes=arrayfun(@(x) number(x.node),source);
        out.weights=arrayfun(@(x) number(x.weight),source);
    end
    function value=number(x)
        if isnumeric(x), value=double(x);
        else, value=str2double(regexprep(x,'`.*$',''));
        end
    end
end
