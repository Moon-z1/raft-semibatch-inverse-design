function d=control_dimension(c)
    d=2*(c.segments-1);
    if isfield(c,'optimize_initial') && c.optimize_initial
        d=d+1;
    end
    if isfield(c,'late_tail') && c.late_tail
        d=d+2;
    end
end
