% load comparison_abc_SaHABC;

t = 500:10:output{1}(end,1)-100;
nt = length(t);
objMeanCurve = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output));
    for jj = 1:length(output)
        ind = find(output{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output{jj}(ind,2);
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

objMeanCurve_SaHABC = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_SaHABC));
    for jj = 1:length(output_SaHABC)
        ind = find(output_SaHABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_SaHABC{jj}(ind,2);
    end
    objMeanCurve_SaHABC(kk) = mean(obj_temp);
end


minObj =  min(min(min(Obj),min(Obj_ABC)),min(Obj_SaHABC));
figure;

semilogy(t,(objMeanCurve-minObj),'k-',t,(objMeanCurve_ABC-minObj),'k:',t,(objMeanCurve_SaHABC-minObj),'k--');


legend('AHABC','ABC','SAHABC');
title('Objective Convergence Speed');
xlabel('Number of function evaluations');
ylabel('The mean of max(lnL)-lnL in 30 runs');


