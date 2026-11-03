function apply_report_style(f)
% 统一报告绘图格式；不改坐标数据、坐标范围和百分数口径。
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
        obj.FontSize=10;
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
