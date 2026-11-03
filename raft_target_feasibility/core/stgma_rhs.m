function dy = stgma_rhs(~,y,p)
% Concentration states, mol/L, constant-volume batch prototype.
% 1:2 monomers; 3 initiator; 4 small radical; 5 original RAFT agent;
% 6:11 active moments [Y0A,Y1A,Y2A,Y0B,Y1B,Y2B];
% 12:17 dormant moments; 18:20 dead moments; 21:22 incorporated units;
% 23 cumulative chain births; 24 cumulative combination events (audit states).
% Reaction assembly corresponds to Tables 1,3,4 of the cited paper.
% Table 4 Eq.16 prints Z1b in one reverse cross-exchange gain; use Z1a
% as required by Table 3 Eq.6, Table 4 Eq.18 and chain-length conservation.
% No adduct fragmentation, intermediate cross termination, solvent transfer,
% propagation diffusion or unreferenced legacy coefficients are retained.
assert(strcmp(p.autoacceleration_mode,'off_low_conversion_only'),...
    'High-conversion onset/free-volume implementation has not been validated.');
M=y(1:2);I=y(3);R=y(4);T=y(5);
Y=reshape(y(6:11),3,2);Z=reshape(y(12:17),3,2);
birth=p.ki.*M*R;
prop=M.*(p.kp.'*Y(1,:).');
shift=[Y(1,:);Y(2,:)+Y(1,:);Y(3,:)+2*Y(2,:)+Y(1,:)];
dY=repmat(birth.',3,1)+(shift*p.kp).*M.'-Y.*(p.kp*M).';
% Reversible activation: T+Y_i <-> R+Z_i.
activation=Y.*(p.kaa*T).'-Z.*(p.kfa*R).';
dY=dY-activation;
dZ=activation;
% RAFT exchange conserves end identity and all total polymer moments.
for i=1:2
    for j=1:2
        % For i~=j, active i becomes dormant i and dormant j becomes active j.
        % Loss and gain moments carry their original chain lengths.
        dY(:,i)=dY(:,i)-p.exchange(i,j)*Y(:,i)*Z(1,j);
        dZ(:,i)=dZ(:,i)+p.exchange(i,j)*Y(:,i)*Z(1,j);
        dY(:,j)=dY(:,j)+p.exchange(i,j)*Z(:,j)*Y(1,i);
        dZ(:,j)=dZ(:,j)-p.exchange(i,j)*Z(:,j)*Y(1,i);
    end
end
kt=p.ktc+p.ktd;
dY=dY-Y.*(kt*Y(1,:).').';
dD=zeros(3,1);
for q=0:2
    dD(q+1)=Y(q+1,:)*p.ktd*Y(1,:).';
    for k=0:q
        dD(q+1)=dD(q+1)+.5*nchoosek(q,k)*Y(k+1,:)*p.ktc*Y(q-k+1,:).';
    end
end
dy=zeros(24,1);
dy(1:2)=-birth-prop;
dy(3)=-p.kd*I;
dy(4)=2*p.f*p.kd*I-sum(birth)+sum(activation(1,:));
dy(5)=-sum(activation(1,:));
dy(6:11)=dY(:);dy(12:17)=dZ(:);dy(18:20)=dD;
dy(21:22)=birth+prop;
dy(23)=sum(birth);
dy(24)=.5*Y(1,:)*p.ktc*Y(1,:).';
end
