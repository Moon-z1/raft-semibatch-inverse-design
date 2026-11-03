% 计算活性链、休眠链及死链的矩方程；状态排列与原模型保持一致。
function dy = stgma_rhs(~,y,p)
    % 恒容间歇模型；浓度单位为mol/L。
    % 状态1:2为单体，3为引发剂，4为小分子自由基，5为初始RAFT试剂。
    % 状态6:11为两类活性链的零、一、二阶矩：[Y0A,Y1A,Y2A,Y0B,Y1B,Y2B]。
    % 状态12:17为休眠链矩，18:20为死链矩，21:22为已聚合单体量。
    % 状态23记录链生成总量，24记录偶合终止总量，用于守恒检查。
    % 反应项对应所引文献表1、3、4。
    % 文献表4式16的反向交叉交换生成项印为Z1b，此处采用Z1a；
    % 依据为文献表3式6、表4式18及链长守恒。
    % 未加入加合物裂解、中间体交叉终止、溶剂链转移，
    % 也未加入增长扩散修正或无文献依据的旧参数。
    assert(strcmp(p.autoacceleration_mode,'off_low_conversion_only'),...
    'High-conversion onset/free-volume implementation has not been validated.');
    M=y(1:2);
    I=y(3);
    R=y(4);
    T=y(5);
    Y=reshape(y(6:11),3,2);
    Z=reshape(y(12:17),3,2);
    birth=p.ki.*M*R;
    prop=M.*(p.kp.'*Y(1,:).');
    shift=[Y(1,:);Y(2,:)+Y(1,:);Y(3,:)+2*Y(2,:)+Y(1,:)];
    dY=repmat(birth.',3,1)+(shift*p.kp).*M.'-Y.*(p.kp*M).';
    % 可逆活化：T+Y_i <-> R+Z_i。
    activation=Y.*(p.kaa*T).'-Z.*(p.kfa*R).';
    dY=dY-activation;
    dZ=activation;
    % RAFT交换保持链端类型及聚合物总矩守恒。
    for i=1:2
        for j=1:2
            % i不等于j时，活性链i转为休眠链i，休眠链j转为活性链j。
            % 损失项与生成项均保留对应链的原链长。
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
    dy(6:11)=dY(:);
    dy(12:17)=dZ(:);
    dy(18:20)=dD;
    dy(21:22)=birth+prop;
    dy(23)=sum(birth);
    dy(24)=.5*Y(1,:)*p.ktc*Y(1,:).';
end
