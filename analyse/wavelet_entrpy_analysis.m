% Quick analysis to see what the wavelet entropy does after the transition
% (to make sure it reaches a maximum or at least stabilises after the
% transition)


load('exp_after_transition.mat')

entropy_vec = zeros(size(sol,2));
for n=1:size(sol,2)
    entropy_vec(n) = wentropy(sol(:,n),'shannon');
end

plot(entropy_vec)