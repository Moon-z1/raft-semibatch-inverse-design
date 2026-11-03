function x=sample_controls(n,c)
D=control_dimension(c);x=zeros(n,D);i=0;
while i<n
 z=2*rand(1,D)-1;[~,~,ok]=feed_decode(z,c);
 if ok,i=i+1;x(i,:)=z;end
end
end
