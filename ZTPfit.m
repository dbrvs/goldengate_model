
%% 
%Daniel Reeves 08/08/2017

%A simulation of drawing 2 viral samples from N people who
%might have R strains of the virus. This R is assumed to come from a ZTP 
%distribution

clf; clear all;

%parameters the user can/should change
averages=10;
N=489; %the size of the cohort 
N2=18; %number superinfected

lambda=linspace(0.01,0.5,100); %the range of ZTP parameters that lead to average richness, low for most cases

%loop over lambda values
for in=1:length(lambda)

l=lambda(in);
avgZTP(in)=l*exp(l)/(exp(l)-1);
stdZTP(in)=sqrt(avgZTP(in)*(1-avgZTP(in)/exp(l)));

y=poissrnd(l,1e6,1); y=y(y>0); %throw out people who have no types to simulate the ZTP

%Generate 2 random vectors with the same probabilities:
cs=0;
%loop over the number of patients in trial but include multiple averages for the trial
for i=1:N*averages

probs = ones(1,y(i)); %flat distribution
%probs = rand(1,y(i)); %random distributed probs
%a=2; probs = exp(-a*flip(1:y(i))); %exponentially decaying prob of viruses

p = probs./sum(probs); %normalize

trial = mnrnd(2,p,1); %do a measurement

if isempty(trial(trial>1));
    cs=cs+1; %counts number measured with the same type twice
end

end

y_dat=[N-N2,N2]*averages;
y_sim=[N*averages-cs,cs];

MSE(in)=sqrt(mean((y_dat-y_sim).^2));  %compute RMS error
end

clear cs i in l p probs trial y y_sim

%% now make the plots
subplot(121)
plot(lambda,MSE,'o','Color','k')
ylim([-1 max(MSE)+1])
xlabel('model parameter \lambda')
ylabel('model score')

subplot(122)
blam=lambda(find(MSE==min(MSE)));
k=1:3;
ZTPpdf = @(x,lam) poisspdf(x,lam) ./ (1-poisscdf(0,lam));

bar([1 2]-0.2,y_dat/N/averages,0.4,'FaceColor',[0 .5 .5])
for ik=1:length(blam)
    hold on
    bar(k+0.2,ZTPpdf(k,blam(ik)),0.4)
    hold off
end
ylim([-0.05,1.1]) 

xlim([.5 3.5])
set(gca,'xTick',[1,2,3])
xlabel('number of strains (n)')
ylabel('probability of n strains')
legend({'data','best-fit model'},'Fontsize',8)
cs=findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',10); end; end

print -dpdf ZTPfit.pdf

avgR=mean(avgZTP(find(MSE==min(MSE))));
stdR=std(avgZTP(find(MSE==min(MSE))));

toc;

disp(['Best \lambda = ' num2str(mean(blam))...
    ' corresponding <R> = ' ...
    num2str(avgR) ' +/- ' num2str(stdR) ])

