function newx=embed_segments(x,oldc,newc)
    assert(newc.segments==2*oldc.segments,'This nesting test doubles segments');
    [~,~,~,w]=feed_decode(x,oldc);
    wn=repelem(w/2,1,2);
    z=log(wn(:,1:end-1)./wn(:,end))/newc.logit_bound;
    newx=[z(1,:),z(2,:)];
    assert(all(abs(newx)<=1+1e-8));
end
