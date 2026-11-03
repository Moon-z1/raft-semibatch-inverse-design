% 按进料时间段积分；切换流量时延续上一段的终点状态。
function o=simulate_reactor(c,rates)
    % 在流量切换处结束本段积分，并以终点状态开始下一段。
    y0=zeros(24,1);
    y0(1:2)=c.n0;
    y0(3)=c.I0;
    y0(5)=c.CTA0;
    [~,c.Vf0]=reactor_properties(y0,c);
    if isempty(rates)
        rates=zeros(2,1);
        edges=[0 c.total_time];
    else
        if isfield(c,'feed_edges')
            edges=[c.feed_edges,c.total_time];
        else
            edges=[linspace(0,c.feed_time,size(rates,2)+1),c.total_time];
        end
        rates=[rates,zeros(2,1)];
    end
    opts=odeset('RelTol',c.RelTol,'AbsTol',c.AbsTol,'MaxStep',c.max_step);
    tt=[];
    yy=[];
    ff=[];
    fed=[];
    current=y0;
    sumfeed=zeros(2,1);
    for j=1:size(rates,2)
        sol=ode15s(@(t,y)reactor_rhs(t,y,c,rates(:,j)),edges(j:j+1),current,opts);
        if sol.x(end)<edges(j+1)-1e-6
            error('Incomplete integration');
        end
        grid=linspace(edges(j),edges(j+1),max(30,ceil((edges(j+1)-edges(j))/c.output_dt)+1));
        if j==1
            grid=unique([grid,logspace(-8,log10(min(c.output_dt,edges(j+1))),80)]);
        end
        if j>1
            grid=grid(2:end);
        end
        vals=deval(sol,grid).';
        tt=[tt;grid.'];
        yy=[yy;vals];
        ff=[ff;repmat(rates(:,j).',numel(grid),1)];
        
        fed=[fed;sumfeed.'+(grid.'-edges(j))*rates(:,j).'];
        
        sumfeed=sumfeed+diff(edges(j:j+1))*rates(:,j);
        current=sol.y(:,end);
    end
    V=zeros(size(tt));
    vf=V;
    dc=V;
    instant=V;
    for k=1:numel(tt)
        [dy,V(k),vf(k),dc(k)]=reactor_rhs(tt(k),yy(k,:).',c,ff(k,:).');
        instant(k)=dy(21)/max(sum(dy(21:22)),realmin);
    end
    B=yy(:,21:22);
    total=c.n0+sumfeed;
    nt=sum(total);
    Xmol=sum(B,2)/nt;
    Xmass=B*c.p.MW/(total.'*c.p.MW);
    Fav=B(:,1)./max(sum(B,2),realmin);
    mom=[sum(yy(:,[6 9 12 15 18]),2),sum(yy(:,[7 10 13 16 19]),2),sum(yy(:,[8 11 14 17 20]),2)];
    DPn=mom(:,2)./max(mom(:,1),realmin);
    DPw=mom(:,3)./max(mom(:,2),realmin);
    % Mn采用矩方程中的重复单元质量，不含端基质量。
    Mn=B*c.p.MW./max(mom(:,1),realmin);
    Mw_approx=DPw.*(Fav*c.p.MW(1)+(1-Fav)*c.p.MW(2));
    checks.component=max(abs(yy(:,1:2)+B-c.n0.'-fed),[],'all')/nt;
    checks.raft=max(abs(yy(:,5)+yy(:,12)+yy(:,15)-c.CTA0))/c.CTA0;
    checks.units=max(abs(mom(:,2)-sum(B,2)))/nt;
    checks.chains=max(abs(mom(:,1)-yy(:,23)+yy(:,24)))/max(mom(:,1));
    checks.minimum_state=min(yy,[],'all');
    checks.volume_fraction=max(abs(sum(fractions(yy,c),2)-1));
    assert(max([checks.component checks.raft checks.units checks.chains])<1e-5 && checks.minimum_state>-1e-7,'Conservation/positivity failed');
    assert(all(DPw(sum(B,2)>1e-7)>=DPn(sum(B,2)>1e-7)-1e-5),'Invalid moments');
    o=struct('t',tt,'y',yy,'feed',ff,'fed',fed,'V',V,'Vf',vf,'dc',dc,'Xmol',Xmol,'Xmass',Xmass,...
    'Finstant',instant,'Fav',Fav,'DPn',DPn,'DPw',DPw,'Mn_repeat',Mn,'Mw_approx',Mw_approx,'checks',checks,'c',c);
end
function phi=fractions(y,c)
    phi=zeros(size(y,1),4);
    for k=1:size(y,1)
        [~,~,phi(k,:)]=reactor_properties(y(k,:).',c);
    end
end
