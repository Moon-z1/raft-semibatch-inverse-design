function apply_report_style(f)
% 统一报告绘图格式；不改坐标数据、坐标范围和百分数口径。
f.Units='pixels';
pos=f.Position;
f.Position=[pos(1:2),1100,1100*pos(4)/pos(3)];
axesList=findall(f,'Type','axes');
for i=1:numel(axesList)
    ax=axesList(i);
    if strcmp(ax.Visible,'off')
        continue
    end
    ax.Box='on';
    ax.LineWidth=0.8;
    ax.XColor='k';
    ax.YColor='k';
    ax.FontName='Times New Roman';
    ax.FontSize=11;
    ax.TickDir='in';
    ax.TickLabelInterpreter='tex';
    if strcmp(ax.XTickLabelMode,'manual')
        ax.XTickLabel=format_text(ax.XTickLabel);
    end
    if strcmp(ax.YTickLabelMode,'manual')
        ax.YTickLabel=format_text(ax.YTickLabel);
    end
end
objects=findall(f);
for i=1:numel(objects)
    obj=objects(i);
    if isprop(obj,'FontName')
        obj.FontName='Times New Roman';
    end
    if isprop(obj,'Interpreter')
        obj.Interpreter='tex';
    end
    if isprop(obj,'String') && ~isempty(obj.String)
        obj.String=format_text(obj.String);
    end
    if strcmp(obj.Type,'legend')
        obj.Box='off';
        obj.FontSize=11;
    elseif strcmp(obj.Type,'text')
        obj.FontSize=12;
    end
end
end

function out=format_text(value)
% 中文用宋体，字母、数字及单位用Times New Roman。
items=string(value);
for i=1:numel(items)
    s=char(items(i));
    s=regexprep(s,'\\fontname\{[^}]*\}','');
    s=strrep(s,'Candidate selection and other','选点等开销');
    s=strrep(s,'Mean online time per run / s','平均在线耗时 / s');
    s=strrep(s,'RMSE reduction relative to scalar GP (%)','相对标量GP的误差降幅 / %');
    s=strrep(s,'(d) New eight targets: equal 30-call budget','(d) 八条新增目标：追加30次评价');
    s=strrep(s,'(c) Per-target improvement','(c) 各目标的误差改善');
    s=strrep(s,'(b) All three repeats retained','(b) 三次重复的最大偏差');
    s=strrep(s,'(a) Target means and ranges across repeats','(a) 三次重复的均值与范围');
    s=strrep(s,'Independent hidden target','独立可达目标编号');
    s=strrep(s,'Maximum absolute error','最大绝对偏差');
    s=strrep(s,'Tracking RMSE','组成RMSE');
    s=strrep(s,'Forward ODE','动力学计算');
    s=strrep(s,'GP fitting','GP拟合');
    s=strrep(s,'Scalar GP','标量GP');
    s=strrep(s,'Trajectory GP','轨迹GP');
    s=strrep(s,'RMSE threshold','RMSE达标线');
    s=strrep(s,'Fixed initial recipe','固定初始配方');
    s=strrep(s,'Mayo-Lewis initial recipe','反算初始配方');
    s=strrep(s,'(c) Gradual gradient','(c) 大幅缓变');
    s=strrep(s,'(b) Steep gradient','(b) 大幅陡变');
    s=strrep(s,'(a) Constant','(a) 恒定组成');
    s=strrep(s,'Instantaneous St fraction','瞬时St摩尔分数');
    s=strrep(s,'Polymerization progress','聚合进度');
    s=strrep(s,'Four stages','四段');
    s=strrep(s,'Ten stages','十段');
    s=strrep(s,'Target','目标');
    s=strrep(s,'Time / h','时间 / h');
    s=strrep(s,'Feed rate','进料流量');
    s=strrep(s,'体积流量','进料流量');
    s=strrep(s,'t / h','时间 / h');
    s=strrep(s,'kappa','κ');
    s=strrep(s,'聚合进度 ξ','聚合进度');
    s=strrep(s,'Spearman排序相关系数','斯皮尔曼秩相关系数');
    s=strrep(s,'原GP-EI','GP-EI');
    s=strrep(s,'均值贪心','GP均值');
    s=strrep(s,'mL s^{-1}','mL\cdot s^{-1}');
    s=strrep(s,'mL·s^{-1}','mL\cdot s^{-1}');
    s=strrep(s,'(s)','/ s');
    % 保留原图中10的幂及希腊字母等TeX表示。
    outText='';
    currentChinese=false;
    for j=1:numel(s)
        code=double(s(j));
        isChinese=(code>=hex2dec('4E00') && code<=hex2dec('9FFF')) ...
            || (code>=hex2dec('3000') && code<=hex2dec('303F')) ...
            || (code>=hex2dec('FF00') && code<=hex2dec('FFEF'));
        if j==1 || isChinese~=currentChinese
            if isChinese
                outText=[outText '\fontname{SimSun}'];
            else
                outText=[outText '\fontname{Times New Roman}'];
            end
        end
        outText=[outText s(j)];
        currentChinese=isChinese;
    end
    items(i)=string(outText);
end
if iscell(value)
    out=cellstr(items);
elseif isstring(value)
    out=items;
else
    out=char(items);
end
end
