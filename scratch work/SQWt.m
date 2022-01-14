A = [15 15 15 15 15]; B = [5 5 5 5 5 5 5 5 5 5]; AB = horzcat(A,B);
ABt = AB;
while length(ABt) < 435
    ABt = [ABt AB];
end
AB = ABt;
for i = 1:length(AB)
     ABt = [ABt;AB];
end
trueSample = ABt'; figure; mesh(trueSample);
s = trueSample(:,randi(425));
h = (linspace(0,length(s),length(s)))'; figure; plot(h,s); xlim([0 length(s)]);

