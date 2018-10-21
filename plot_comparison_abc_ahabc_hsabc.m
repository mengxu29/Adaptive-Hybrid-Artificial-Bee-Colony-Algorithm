load comparison_abc_ahabc_hsabc_18;

t = 500:10:output_AHABC{1}(end,1)-100;
nt = length(t);
objMeanCurve = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_AHABC));
    for jj = 1:length(output_AHABC)
        ind = find(output_AHABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_AHABC{jj}(ind,2);
    end
    objMeanCurve(kk) = mean(obj_temp);
end

objMeanCurve_ABC = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_ABC));
    for jj = 1:length(output_ABC)
        ind = find(output_ABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_ABC{jj}(ind,2);
    end
    objMeanCurve_ABC(kk) = mean(obj_temp);
end

objMeanCurve_hsabc = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_HSABC));
    for jj = 1:length(output_HSABC)
        ind = find(output_HSABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_HSABC{jj}(ind,2);
    end
    objMeanCurve_hsabc(kk) = mean(obj_temp);
end


minObj =  min(min(min(Obj_AHABC),min(Obj_ABC)),min(Obj_HSABC));
figure;

semilogy(t,(objMeanCurve-minObj),'k-',t,(objMeanCurve_ABC-minObj),'r-',t,(objMeanCurve_hsabc-minObj),'k--');


legend('AHABC','ABC','HSABCA');
title('Objective Convergence speed');
xlabel('Number of Function Evaluations');
ylabel('The mean of max(lnL)-lnL in 30 runs');


