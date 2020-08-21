%% load recorded data
time        = 0:23; %[hour] 24-hour clock
temperature = [36 35.8 36.2 35.5 35.8 35.3 35.5 34.2 35.7 36.8 36.1 35.9 ...
               34.1 37.6 37.7 35 34.2 35.5 36.8 35.1 34.9 35.6 36.5 37.1]; %[degree Celsius] day1

%% section "A plot and a table"
figure;
h1 = plot(time,temperature,'bo-'); hold on;
h2 = plot([0 23],38*[1 1],'r--', 'displayname','high temperature');
xlabel('Time [hour]')
ylabel('Temperature [^oC]')
grid on;
ylim([34 39]);
legend(h2, 'location','northeast');

%% section "Maximum, minimum and range"
max(temperature)
min(temperature)
range(temperature)

%% section "Mean, median, mode"
mean(temperature)
sort(temperature)
median(temperature)
mode(temperature)

%% section "Mean, median, mode: visualization"
figure;
h1 = plot(time,temperature,'bo-', 'displayname',''); hold on;
h2 = plot([0 23],mean(temperature)*[1 1],'k-', 'linewidth',2, 'displayname','mean');
h3 = plot([0 23],median(temperature)*[1 1],'r-', 'linewidth',2, 'displayname','median');
h4 = plot([0 23],mode(temperature)*[1 1],'g-', 'linewidth',2, 'displayname','mode');
xlabel('Time [hour]');
ylabel('Temperature [^oC]');
grid on;
ylim([34 39]);
legend([h2 h3 h4]); 


%% section "Boxplot"
figure;
h = boxplot(temperature);
set(h,{'linew'},{2}); 
grid on;
ylabel('Temperature [^oC]');
ylim([34 39]);
xlabel('Day');

%% section "Variance, standard deviation"
figure;
h1 = plot(time,temperature,'bo-', 'displayname','data x'); hold on;
h2 = plot([0 23],mean(temperature)*[1 1],'k-', 'linewidth',2, 'displayname',['mean \mu=',num2str(mean(temperature),"%.1f")]);
h3 = plot([0 23],mean(temperature)+std(temperature)*[1 1],'k--', 'linewidth',2, 'displayname',['\mu+\sigma=',num2str(mean(temperature)+std(temperature),"%.1f")]);
h4 = plot([0 23],mean(temperature)-std(temperature)*[1 1],'--', 'color',.7*[1 1 1], 'linewidth',2, 'displayname',['\mu-\sigma=',num2str(mean(temperature)-std(temperature),"%.1f")]);
xlabel('Time [hour]');
ylabel('Temperature [^oC]');
grid on;
ylim([34 39]);
legend([h1 h2 h3 h4]);

%% section "Variance, standard deviation: usage"
x = temperature;
mu = mean(temperature);
y = mu + .5*(x-mu); %[degree Celsius]

figure;
h1 = plot(time,y,'bo-', 'displayname','data y'); hold on;
h2 = plot([0 23],mean(y)*[1 1],'k-', 'linewidth',2, 'displayname',['mean \mu=',num2str(mean(y),"%.1f")]);
h3 = plot([0 23],mean(y)+std(y)*[1 1],'k--', 'linewidth',2, 'displayname',['\mu+\sigma_y=',num2str(mean(y)+std(y),"%.1f")]);
h4 = plot([0 23],mean(y)-std(y)*[1 1],'--', 'color',.7*[1 1 1], 'linewidth',2, 'displayname',['\mu-\sigma_y=',num2str(mean(y)-std(y),"%.1f")]);
xlabel('Time [hour]');
ylabel('Temperature [^oC]');
grid on;
ylim([34 39]);
legend([h1 h2 h3 h4]);

%% section "The histogram"
x = temperature;

figure;
h1 = plot(time,temperature,'bo-', 'displayname','data x'); hold on;
xlabel('Time [hour]');
ylabel('Temperature [^oC]');
grid on;
ylim([34 39]);
legend(h1);

f=figure;
set(f,'position',[609 342 439 420])
nbins = 6; %choose
hist(x,nbins);
xlim([34 39]);
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
xlabel('Temperature [^oC]')
ylabel('Number of points')
grid on;
legend(['nbins=',num2str(nbins,"%d")], 'location','northwest');

%% section "Probability distribution (PDF)"
%choose
x = linspace(34,39,1e2);
mu = min(x) + range(x)/2;
sigma = 1;

%consequence
y = exp(-.5*power((x-mu)/sigma,2))/(sigma*sqrt(2*pi));

f = figure;
f.Position=[f.Position(1:2) 238*ones(1,2)];
plot(x,y); grid on;
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
xlabel('Temperature x [^oC]');
ylabel("P(x)");
title("PDF: Normal distribution");
axis tight; 

%% slide "Normal (Gaussian) distribution"
figure;
mu=0;
sigma=1;
plot(x,normpdf(x,mu,sigma), 'linewidth',2, 'displayname',['Normal (Gaussian),\newline{} \mu=',num2str(mu,"%.1f"),' (mean),\newline{} \sigma=',num2str(sigma,"%.1f"), ' (std dev)']); hold on;
grid on;
xlabel('x'); ylabel('P(x)');
%title("Normal (Gaussian) distribution");
legend('show', 'location','northwest');

%% slide "Student's t-distribution"
figure;
nu_dof = 5; %choose
plot(x,tpdf(x,nu_dof), 'linewidth',2, 'displayname',['Student`s t, \nu=',num2str(nu_dof,"%d"),' (dof)']); hold on;
grid on;
xlabel('x'); ylabel('P(x)');
%title("Probability distributions (PDFs)");
legend('show');
ylim([0 .42]);


figure;
nu_dof=1; h1 = plot(x,tpdf(x,nu_dof), 'linewidth',2, 'displayname',['Student`s t, \nu=',num2str(nu_dof,"%d"),' (dof)']); hold on;
nu_dof=3; h2 = plot(x,tpdf(x,nu_dof), 'linewidth',2, 'displayname',['Student`s t, \nu=',num2str(nu_dof,"%d"),' (dof)']); hold on;
nu_dof=5; h3 = plot(x,tpdf(x,nu_dof), 'linewidth',2, 'displayname',['Student`s t, \nu=',num2str(nu_dof,"%d"),' (dof)']); hold on;
          h4 = plot(x,normpdf(x,0.,1.),'k-', 'linewidth',2, 'displayname',['Normal (Gaussian),\newline{} \mu=',num2str(mu,"%.1f"),' (mean),\newline{} \sigma=',num2str(sigma,"%.1f"), ' (std dev)']); hold on;
grid on;
xlabel('x'); ylabel('P(x)');
%title("Probability distributions (PDFs)");
legend([h1 h2 h3 h4]);
ylim([0 .42]);



%% section "\chi^2 distribution"
x = linspace(-5,5,1e2); %choose

mu = 0; %mean
sigma = 1;
y_npdf = normpdf(x,mu,sigma);

nu_dof = 3; %choose
y_tpdf = tpdf(x,nu_dof);

k_dof = 2.1;
y_chiSqpdf = chi2pdf(x,k_dof);

%%%% compared study of multiple distributions
figure;
plot(x,y_npdf, 'linewidth',2, 'displayname',['Normal (Gaussian),\newline{} \mu=',num2str(mu,"%.1f"),' (mean),\newline{} \sigma=',num2str(sigma,"%.1f"), ' (std dev)']); hold on;
plot(x,y_tpdf, 'linewidth',2, 'displayname',['Student t, \nu=',num2str(nu_dof,"%d"),' (dof)']); hold on;
plot(x,y_chiSqpdf, 'linewidth',2, 'displayname',['\chi^2, k=',num2str(k_dof,"%.2f"),' (dof)']); hold on;
grid on;
xlabel('x'); ylabel('P(x)');
title("Probability distributions (PDFs)");
legend('show');


%% section "Correlation"
x1 = temperature; %[degree Celsius]

rng(22); %fix random seed for repeatability purpose 
x2 = x1 + randn(size(x1)); %[degree Celsius]

figure;
h1 = plot(time,x1,'bo-', 'linewidth',2, 'displayname','day1: data x_1'); hold on;
h2 = plot(time,x2,'mo-', 'linewidth',2, 'displayname','day2: data x_2');
xlabel('Time [hour]');
ylabel('Temperature [^oC]');
grid on;
ylim([34 39]);
legend('show', 'location','northwest');

figure;
plot(x1(1:end),x2(1:end),'bo', 'linewidth',2, 'displayname','(x_1,x_2) exp data'); grid on; 
xlabel('x_1 [^oC]');
ylabel('x_2 [^oC]');
axis equal; 
xlim([34 39]);
ylim([34 39]);

%% slide "Simple linear regression"
%WIP

%% slide "Time series"
%%%% histogram corresponding to gamma=2
f=figure;
set(f,'position',[609 342 439 420])
nbins = 6; %choose
hist(x2,nbins);
xlim([34 39]);
set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
xlabel('Temperature [^oC]')
ylabel('Number of points')
grid on;
legend(['nbins=',num2str(nbins,"%d")], 'location','northwest');

%%%% histogram corresponding to time=2:00
gamma_max = 30;
x2D = [x1; ...
       x2; ...
       ones(gamma_max,1).*x + randn(gamma_max,length(x1))] ; %multivariate time series 

f=figure;
set(f,'position',[609 342 439 420])
nbins = 6; %choose
k = 3; %choose
hist(x2D(:,k),nbins);
xlim([34 39]);
%set(gca,'xdir','reverse'); %alternatively use GUI: right click on the axes of the figure> Open Property Inspector > Rulers > XDir > reverse
xlabel('Temperature [^oC]')
ylabel('Number of points')
grid on;
legend(['nbins=',num2str(nbins,"%d")], 'location','northwest');

%%
