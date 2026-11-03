function x=neutral_control(c)
    x=zeros(1,control_dimension(c));
    if isfield(c,'late_tail') && c.late_tail
        d=2*(c.segments-1);
        x(d+(1:2))=-1;
    end
end
