%% Experimentally observe how Chi distribution is/gets formed

%choose
nbins = 2e2;

% nu = v.a.i.i.d. with i.d. standard normal
nu = randn(1e5,1);

figure;
hist(nu,nbins); grid on; 
ylabel('Frequency i.e. nr of occurrences');
xlabel('\nu');
title('Histogram of random variables \nu drawn from a\newline{}parent population that is standard normal');


%% Chi: to be compared with https://en.wikipedia.org/wiki/Chi-square_distribution 
%choose
k = 6;

%ini
q = nan(2e5,1);
for i=1:length(q)
   q(i)= sum(power(randn(k,1),2) ); 
end

figure;
h = histogram(q,nbins);
ylabel('Frequency i.e. nr of occurrences');
xlabel('q');
title(['Histogram of derived random variable q=X_1^2+X_2^2+..+X_',num2str(k,"%d"),'^2,\newline{} where X_1 to X_',num2str(k,"%d"),' form a sample drawn from a\newline{} parent population that is standard normal\newline{}(the distribution of q converges to \chi^2(',num2str(k,"%d"),') on the long-run']);
grid on;

deltax = (max(q)-min(q))/nbins ; %=h.BinWidth
integral_val = sum(h.Values)*h.BinWidth;

figure; %plot(h.Values)
plot(h.BinEdges(1:end-1)+h.BinWidth/2, h.Values/integral_val) ;
grid on;
title(['Normalized histogram of q']);
ylabel('Normalized frequency (or density)');
xlabel('q');
xlim([0 8]);

