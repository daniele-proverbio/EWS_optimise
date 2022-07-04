%% Test significance

function param_ok = testsignificance(matrix,val2)

p_value = zeros(size(matrix,1)-1,1);  % comparison with no change at all
for n=2:size(matrix,1)-1
            [~, pv] = ttest(matrix(size(matrix,1)-n,:) , matrix(size(matrix,1)-1,:));
            p_value(n) = abs(pv);
end

p_value_mean = movmean(p_value,10);  % a bit of smoothing to suppress fluctuations due to finite sizes
position_ok = find(p_value_mean < 0.05);  % when is p-value < 0.05 (typical significance)
check_position = diff(position_ok);
check = 1;
param_ok = val2(position_ok(1));

% find the first value when p-value <0.05 stably (i.e., it is followed by
% other ~10 values <0.05) -> to avoid picking up fluctuations but instead
% being robust
for n =1:size(check_position,1)-10
    for m=0:9
        if check_position(n+m) ~= 1
            check = 0;
            param_ok = val2(position_ok(n+1));
        else
            check = 1;
        end
    end
end

% % % Example figure
% figure()
% hold on
% plot(val2(1:end-1),p_value)
% plot(val2(1:end-1),p_value_mean)
% yline(0.05)
% xline(param_ok)
% legend({'p-value','Moving avg','p_v = 0.05'},FontSize=14,Location='northwest')


end