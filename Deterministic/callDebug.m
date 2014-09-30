t=linspace(0,5000,10000);

P00 = 1;
P01 = 0;
P10 = 0;
P11 = 0;
P12 = 0;
P02 = 0;
marRA = 0;
Auf = 0;
Ruf = 0;
A = 0;
R = 0;
R2 = 0;

MAR0 = [ P00,P01,marRA,R2];

salCon = 0.;



lr= 0.02; %min^-1 Nicoloff 2006. Increased


k_r = 1.;
kr = 0.002;
lmar = 0.02; %min^-1  arep.med.harvard.edu/rna_decay ln(2)/24.1 min (Bundschuh 2003)

br=  100;  % min^-1 Martin 2004
a01 = 0.002;
a00 = 0.002;
tic



% Compute all the corresponding values of N
[T,MARcomputed]=ode23s(@marModelSimpleDebug,t,MAR0,[],la, lr, k_a, ka, k_r, kr, a00, a01, a10, lmar, br, ba);
toc

%MARcomputed
% %% Stochastic
% %
%
% fileToRead1 = 'data.dat';
% rawData1 = importdata(fileToRead1);
% [~,name] = fileparts(fileToRead1);
% newData1.(genvarname(name)) = rawData1;
%
% data = newData1.(genvarname(name));
%
%
% time = data(:,1); % timeR
%
%
% y = data(:,11); % A
% x = data(:,13); % R2
%
%
%
% hold on
% plot(time,y)
% plot(time,x,'red')


%% Deterministic
hold on
plot(T,MARcomputed(:,3),'LineWidth',2)
plot(T,MARcomputed(:,4),'red','LineWidth',2)
% legend('A','R2')
format shortE
MARcomputed(end,:)
