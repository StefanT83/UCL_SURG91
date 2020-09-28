%% slide "PDF of the parent population"
%choose
sigma = 4;
mu = 110;
N = 1e5;
X = mu + sigma*randn(N,1); % X = X(p), where the persons p=[p1 p2 p3 ...];  the measures (i.e. observations) associated to each person: time to complete the track [seconds]

%%%%% plot the histogram: Frequency
nbins = 100; %choose
[nrElemPerBin,posBinCenter] = hist(X,nbins) ;

figure;
h1 = bar(posBinCenter,nrElemPerBin,'c'); grid on;
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
%title('');
xlabel('X [sec]');
ylabel('Frequency (Number of points)');
xlim(mu+4*sigma*[-1 1]);
legend(h1,['histogram of the parent population\newline{}nbins=',num2str(nbins,"%d")] ,'location','northwest');
ylim([0 8000]); 

%%%%% plot the histogram: Relative frequency i.e. a lumped approximation of the discrete distribution of the parent population   
figure;
h1 = bar(posBinCenter,nrElemPerBin/N,'c'); hold on; grid on; %note that N == sum(nrElemPerBin)  
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
%title('');
xlabel('X [sec]');
ylabel('Relative frequency [-]');
xlim(mu+4*sigma*[-1 1]);
legend(h1,['a lumped approximation of the\newline{}discrete distribution of the parent population\newline{}nbins=',num2str(nbins,"%d")], 'location','northwest'); 
ylim([0 .088]);

%%%%% plot the histogram: Scaled relative frequency
area_under_hist = diff(posBinCenter(1:2))*N; %integral under the histogram; note that N == sum(nrElemPerBin)     

figure;
h1 = bar(posBinCenter,nrElemPerBin/area_under_hist,'c'); hold on; grid on;
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
%title('');
xlabel('X [sec]');
ylabel('(Relative frequency)/\Delta{}X');
xlim(mu+4*sigma*[-1 1]);
ylim([0 .12]);
legend(h1,['histogram of the scaled relative frequency\newline{}nbins=',num2str(nbins,"%d")], 'location','northwest'); 


%%%%% plot the histogram: Scaled relative frequency + continuous pdf of the conceptual population 
X_axis = linspace(min(X),max(X),1e2);
y_normPdf = exp(-.5*power((X_axis-mu)/sigma,2))/(sigma*sqrt(2*pi));

figure;
h1 = bar(posBinCenter,nrElemPerBin/area_under_hist,'c'); hold on;
h2 = plot(X_axis,y_normPdf,'b', 'linewidth',4); grid on;
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
title('PDF');
xlabel('X [sec]');
ylabel('P(X)');
xlim(mu+4*sigma*[-1 1]);
ylim([0 .12]);
legend([h1 h2],{['histogram of the scaled relative frequency\newline{}nbins=',num2str(nbins,"%d")], 'continuous PDF of the conceptual parent population'}, 'location','northwest'); 


%% slide "PDF of the sample mean"

%choose sample size
n = 100;

%%% instead of listing all possible samples, here we are only going to consider a few samples formed "at random" from the parent population     
%choose how many samples to generate at random
nrSamples = 1e4;

%preallocate memory
barX = nan(1,nrSamples)

for idx=1:nrSamples
    sample = X(randi(N,1,n)) ; 
    barX(idx) = sum(sample)/n ; % barX = barX(sample)
end

%%%%% plot the histogram
nbins = 100; %choose
[nrElemPerBin,posBinCenter] = hist(barX,nbins) ;

figure;
h1 = bar(posBinCenter,nrElemPerBin/sum(nrElemPerBin),'c'); hold on; grid on;
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
%title('');
xlabel('$$\bar{X}$$ [sec]', 'interpreter','latex');
ylabel('Relative frequency');
xlim(mu+4*sigma*[-1 1]);
legend(h1,['a lumped approximation of the\newline{}discrete distribution of the sample mean\newline{}nbins=',num2str(nbins,"%d"),'; n=',num2str(n,"%d")], 'location','northwest'); 
ylim([0 .04]);

