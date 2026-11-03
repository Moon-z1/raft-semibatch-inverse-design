function [x,c]=apply_design(x,c)
    x=x(:).';
    if isfield(c,'optimize_initial') && c.optimize_initial
        total=c.n0+c.dose;
        f=.70+.15*x(end);
        x=x(1:end-1);
        vm=c.p.MW./(1000*c.p.rho_monomer_g_mL);
        frac=[f;1-f];
        c.n0=frac*c.V0/(frac.'*vm);
        c.dose=total-c.n0;
        assert(all(c.dose>0),'Initial split exceeds total supply');
        c.initial_fraction=f;
        c.optimize_initial=false;
    end
end
