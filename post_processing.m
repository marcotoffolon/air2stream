%-----------------------------------------------------------------------
%               POST-PROCESSING OF AIR2STREAM OUTPUTs
%-----------------------------------------------------------------------

close all
clear all
clc

%------------------------> main informations:

nome_fiume = 'Switzerland';         % name of folder
stazione_a = 'KOP';                 % name/ID of the air station
stazione_w = '2029';                % name/ID of the water station
time_resolution = '1d';             % time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
numero_parametri = '8';             % version: 3,4,5,7,8 parameters                             %
obiettivo = 'RMS';                  % objective function: KGE, NSE, RMS
mod_num = 'CRN';                    % mod_num :   RK4 , EUL , RK2 , CRN
run_mode = 'PSO';                   % optimization algorithm: PSO or GLUE

flag_obj = 'NSE';                   % option for dottyplots: NSE or NSEmod (NSE in relationship with the mean year)

%------------------------
model = 'air2stream';

ERRORI_cal = NaN(28,1);
ERRORI_val = NaN(28,1);

mkdir([nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri '_' obiettivo]);
mkdir([nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution]);
mkdir([nome_fiume '/errors']);

cartella=[ nome_fiume '/' model '/output_' mod_num '_' numero_parametri '/'];
dati_obiettivo = load([cartella '1_' run_mode '_' obiettivo '_' stazione_a '_' stazione_w '_c_' time_resolution '.out']);
dati_cal = load([cartella '2_' run_mode '_' obiettivo '_' stazione_a '_' stazione_w '_cc_' time_resolution '.out']);
dati_val = load([cartella '3_' run_mode '_' obiettivo '_' stazione_a '_' stazione_w '_cv_' time_resolution '.out']);
primo_val_c = find(dati_cal(:,1) ~= -999);
primo_val_c = primo_val_c(1);
dati_cal = dati_cal(primo_val_c:end,:);
primo_val_v = find(dati_val(:,1) ~= -999);
primo_val_v = primo_val_v(1);
dati_val = dati_val(primo_val_v:end,:);

%---------------------------> read file 0_PSO

fid=fopen([cartella '0_' run_mode '_' obiettivo '_' stazione_a '_' stazione_w '_c_' time_resolution '.out']);
prec='double';                          % precision of file 0_PSO (double)
data=fread(fid,prec);
data=reshape(data,9,length(data)/9)';   % matrix of file 0_PSO
toll = 2;
fclose(fid);
if strcmp(obiettivo,'RMS')
    data(:,end)=-data(:,end);
end

I=find(data(:,end)<=toll);
B=data(I,:);            % B is the name of the file saved in file 0_PSO

a1 = B(:,1);            % values of parameters at each iteration
a2 = B(:,2);
a3 = B(:,3);
a4 = B(:,4);
a5 = B(:,5);
a6 = B(:,6);
a7 = B(:,7);
a8 = B(:,8);
obj_RMS = B(:,9);      % values of efficiency at each iteration


%------------------------> efficiency in calibration and validatio period

obiettivo_cal = dati_obiettivo(9);
obiettivo_val = dati_obiettivo(10);

%------------------------------------------------------------------

dati_cal(dati_cal == -999) = NaN;
dati_val(dati_val == -999) = NaN;

time_calib = datenum(dati_cal(:,1),dati_cal(:,2),dati_cal(:,3));
Ta_calib = dati_cal(:,4);
Tw_oss_calib = dati_cal(:,5);
Tw_sim_d_calib = dati_cal(:,6);
Tw_oss_aggr_calib = dati_cal(:,7);
Tw_sim_aggr_calib = dati_cal(:,8);
Q_calib = dati_cal(:,9);

time_valid = datenum(dati_val(:,1),dati_val(:,2),dati_val(:,3));
Ta_valid = dati_val(:,4);
Tw_oss_valid = dati_val(:,5);
Tw_sim_d_valid = dati_val(:,6);
Tw_oss_aggr_valid = dati_val(:,7);
Tw_sim_aggr_valid = dati_val(:,8);
Q_valid = dati_val(:,9);

%----------------------> mean year of calibration

nd=length(Ta_calib);
for i=1:365
    ii=floor(i:365.25:nd);
    ppa=Ta_calib(ii);         ppw_oss=Tw_oss_calib(ii);   ppw_sim=Tw_sim_d_calib(ii);  ppq=Q_calib(ii);
    ppa(isnan(ppa))=[];    ppw_oss(isnan(ppw_oss))=[];   ppw_sim(isnan(ppw_sim))=[];  ppq(isnan(ppq))=[];
    TAam_calib(i)=mean(ppa);     TWam_oss_calib(i)=nanmean(ppw_oss);  TWam_sim_calib(i)=nanmean(ppw_sim);     Qam_calib(i)=nanmean(ppq);
end

%----------------------> mean year of validation

nd=length(Ta_valid);
for i=1:365
    ii=floor(i:365.25:nd);
    ppa=Ta_valid(ii);         ppw_oss=Tw_oss_valid(ii);   ppw_sim=Tw_sim_d_valid(ii);  ppq=Q_valid(ii);
    ppa(isnan(ppa))=[];    ppw_oss(isnan(ppw_oss))=[];   ppw_sim(isnan(ppw_sim))=[];  ppq(isnan(ppq))=[];
    TAam_valid(i)=nanmean(ppa);     TWam_oss_valid(i)=nanmean(ppw_oss);  TWam_sim_valid(i)=nanmean(ppw_sim);     Qam_valid(i)=nanmean(ppq);
end


dotty_am = [];
years = year(time_calib(1)):year(time_calib(end));
for i=1:length(years)
    dotty_am = [dotty_am TWam_oss_calib];
    if leapyear(years(i))
        dotty_am = [dotty_am TWam_oss_calib(end)];
    end
    
end

dotty_am = dotty_am';
I = find(isnan(Tw_oss_calib)==0);
numero = length(I);
if strcmp(flag_obj,'NSE')
    aaa = (Tw_oss_calib(I) - mean(Tw_oss_calib(I))).^2;
elseif strcmp(flag_obj,'NSEmod')
    aaa = (Tw_oss_calib(I) - dotty_am(I)).^2;
end
sss = sum(aaa);
obj = 1-(numero*(obj_RMS).^2/sss);

if strcmp(flag_obj,'NSE')
    obj_max = max(obj) + 0.001;
    min_plot = obj_max-0.05;
elseif strcmp(flag_obj,'NSEmod')
    obj_max = max(obj)+0.01;
    min_plot = max(obj)-0.5;
end

%------------------------> DOTTYPLOTS

if strcmp(numero_parametri,'3')==1
    if sum(a4)~=0 | sum(a5)~=0 | sum(a6)~=0 | sum(a7)~=0 | sum(a8)~=0
        disp('errore in calibrazione');
        break
    else
        figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 25 7.5]);
        subplot(1,3,1); plot(a1,obj,'.')
        grid on
        axis([-5 15 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a1','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(1,3,2); plot(a2,obj,'.')
        grid on
        axis([0 1.5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a2','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(1,3,3); plot(a3,obj,'.')
        grid on
        axis([-5 5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a3','FontSize',10);
        ylabel('objective function','FontSize',10);
        name1 = ['Dottyplots ' stazione_a '-' stazione_w ' ' time_resolution ' ' flag_obj ' npar=' num2str(numero_parametri)];
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,name1,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',10);
        
        fpat=[nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution '/'];
        
        set(gcf,'renderer', 'zbuffer')
        fnamDOT=['dottyplots' stazione_a '_' stazione_w '_' time_resolution '_' flag_obj '.png'];
        print('-dpng', fullfile(fpat, fnamDOT));
        
        close all
    end
end

if strcmp(numero_parametri,'4')==1
    if sum(a5)~=0 | sum(a6)~=0 | sum(a7)~=0 | sum(a8)~=0
        disp('errore in calibrazione');
        break
    else
        figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 25 15]);
        subplot(2,3,1); plot(a1,obj,'.')
        grid on
        axis([0 10 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a1','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,2); plot(a2,obj,'.')
        grid on
        axis([0 1.5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a2','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,3); plot(a3,obj,'.')
        grid on
        axis([-5 5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a3','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,4); plot(a4,obj,'.')
        grid on
        axis([-1 1 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a4','FontSize',10);
        ylabel('objective function','FontSize',10);
        name1 = ['Dottyplots ' stazione_a '-' stazione_w ' ' time_resolution ' ' flag_obj ' npar=' num2str(numero_parametri)];
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,name1,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',10);
        
        fpat=[nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution '/'];
        
        set(gcf,'renderer', 'zbuffer')
        fnamDOT=['dottyplots' stazione_a '_' stazione_w '_' time_resolution '_' flag_obj '.png'];
        print('-dpng', fullfile(fpat, fnamDOT));
        close all
    end
end

if strcmp(numero_parametri,'5')==1
    if sum(a4)~=0 | sum(a5)~=0 | sum(a8)~=0
        disp('errore in calibrazione');
        break
    else
        figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 25 15]);
        subplot(2,3,1); plot(a1,obj,'.')
        grid on
        axis([0 10 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a1','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,2); plot(a2,obj,'.')
        grid on
        axis([0 1.5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a2','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,3); plot(a3,obj,'.')
        grid on
        axis([-5 5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a3','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,5); plot(a6,obj,'.')
        grid on
        axis([0 10 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a6','FontSize',10);
        ylabel('objective function','FontSize',10);
        subplot(2,3,6); plot(a7,obj,'.')
        grid on
        axis([0 1 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a7','FontSize',10);
        ylabel('objective function','FontSize',10);
        name1 = ['Dottyplots ' stazione_a '-' stazione_w ' ' time_resolution ' ' flag_obj ' npar=' num2str(numero_parametri)];
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,name1,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',10);
        
        fpat=[nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution '/'];
        
        set(gcf,'renderer', 'zbuffer')
        fnamDOT=['dottyplots' stazione_a '_' stazione_w '_' time_resolution '_' flag_obj '.png'];
        print('-dpng', fullfile(fpat, fnamDOT));
        close all
    end
end

if strcmp(numero_parametri,'7')==1
    if sum(a4)~=0
        disp('errore in calibrazione');
        break
    else
        figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 24 12]);
        subplot(2,4,1); plot(a1,obj,'.')
        grid on
        axis([0 10 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a1','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,2); plot(a2,obj,'.')
        grid on
        axis([0 1.5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a2','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,3); plot(a3,obj,'.')
        grid on
        axis([-5 5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a3','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,5); plot(a5,obj,'.')
        grid on
        axis([0 20 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a5','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,6); plot(a6,obj,'.')
        grid on
        axis([0 10 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a6','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,7); plot(a7,obj,'.')
        grid on
        axis([0 1 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a7','FontSize',8);
        ylabel('objective function','FontSize',8);
        subplot(2,4,8); plot(a8,obj,'.')
        grid on
        axis([-1 5 min_plot obj_max])
        set(gca, 'FontSize', 7);
        xlabel('a8','FontSize',8);
        ylabel('objective function','FontSize',8);
        name1 = ['Dottyplots ' stazione_a '-' stazione_w ' ' time_resolution ' ' flag_obj ' npar=' num2str(numero_parametri)];
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,name1,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',8);
        
        fpat=[nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution '/'];
        
        set(gcf,'renderer', 'zbuffer')
        fnamDOT=['dottyplots' stazione_a '_' stazione_w '_' time_resolution '_' flag_obj '.png'];
        print('-dpng', fullfile(fpat, fnamDOT));
        close all
    end
end

if strcmp(numero_parametri,'8')==1
    figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 24 12]);
    subplot(2,4,1); plot(a1,obj,'.')
    grid on
    axis([0 10 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a1','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,2); plot(a2,obj,'.')
    grid on
    axis([0 1.5 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a2','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,3); plot(a3,obj,'.')
    grid on
    axis([-5 5 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a3','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,4); plot(a4,obj,'.')
    grid on
    axis([-1 1 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a4','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,5); plot(a5,obj,'.')
    grid on
    axis([0 20 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a5','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,6); plot(a6,obj,'.')
    grid on
    axis([0 10 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a6','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,7); plot(a7,obj,'.')
    grid on
    axis([0 1 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a7','FontSize',8);
    ylabel('objective function','FontSize',8);
    subplot(2,4,8); plot(a8,obj,'.')
    grid on
    axis([-1 5 min_plot obj_max])
    set(gca, 'FontSize', 7);
    xlabel('a8','FontSize',8);
    ylabel('objective function','FontSize',8);
    name1 = ['Dottyplots ' stazione_a '-' stazione_w ' ' time_resolution ' ' flag_obj ' npar=' num2str(numero_parametri)];
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,name1,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',8);
    
    fpat=[nome_fiume '/' model '/figure_' run_mode '_' mod_num '_' numero_parametri  '_' obiettivo '/' time_resolution '/'];
    
    set(gcf,'renderer', 'zbuffer')
    fnamDOT=['dottyplots' stazione_a '_' stazione_w '_' time_resolution '_' flag_obj '.png'];
    print('-dpng', fullfile(fpat, fnamDOT));
    close all
end




Tmax = max(max(Ta_calib),max(Tw_oss_calib));
Tmaxx = 1.15*Tmax;
Tmin = min(min(Ta_calib),min(Tw_oss_calib));
Tminn = 1.15*Tmin;
bisettrice = [ Tminn Tmaxx];

Tmax_med = max(max(TAam_calib),max(TWam_oss_calib));
Tmaxx_med = 1.15*Tmax_med;
Tmin_med = min(min(TAam_calib),min(TWam_sim_calib));
Tminn_med = 1.15*Tmin_med;
bisettrice_med = [ Tminn_med Tmaxx_med];


%----------------------------------------------------------------------
%-------------------------> ERRORS in calibration period for daily data
%----------------------------------------------------------------------

RMS_cal = (nanmean( (Tw_oss_calib-Tw_sim_d_calib).^2 ))^0.5;

MAE_cal = nanmean(abs(Tw_oss_calib-Tw_sim_d_calib));

std_obs_cal = nanstd(Tw_oss_calib);
std_mod_cal = nanstd(Tw_sim_d_calib);
mean_obs_cal = nanmean(Tw_oss_calib);
mean_mod_cal = nanmean(Tw_sim_d_calib);
covar_mod_cal = nancov(Tw_oss_calib,Tw_sim_d_calib);
covar_mod_cal = covar_mod_cal(1,2);

alpha_cal = std_mod_cal/std_obs_cal;
beta_cal = mean_mod_cal/mean_obs_cal;
gamma_cal = covar_mod_cal/(std_obs_cal*std_mod_cal);
KGE_cal = 1- sqrt((std_mod_cal/std_obs_cal-1)^2+(mean_mod_cal/mean_obs_cal-1)^2+...
    (covar_mod_cal/(std_mod_cal*std_obs_cal)-1)^2);


sopra_cal = [NaN(length(Tw_oss_calib),1)];
sotto_cal = [NaN(length(Tw_oss_calib),1)];
sopra_cal = nansum((Tw_oss_calib - Tw_sim_d_calib).^2);
sotto_cal = nansum((Tw_oss_calib - mean_obs_cal).^2);
NSE_cal = 1-(sopra_cal/sotto_cal);


% %----> errors at weekly time scale but with calibration at daily time scale

L1w = floor(length(Tw_oss_calib)/7);
TW_oss_1w = [NaN(L1w,1)];
TW_mod_1w = [NaN(L1w,1)];
p=0;
for i = 1:L1w
    TW_oss_1w(i)=nanmean(Tw_oss_calib(p+1:p+7));
    TW_mod_1w(i)=nanmean(Tw_sim_d_calib(p+1:p+7));
    p=p+7;
end

RMS_agg_1w = (nanmean( (TW_oss_1w-TW_mod_1w).^2 ))^0.5;

MAE_agg_1w = nanmean(abs(TW_oss_1w-TW_mod_1w));

std_obs_agg_1w = nanstd(TW_oss_1w);
std_mod_agg_1w = nanstd(TW_mod_1w);
mean_obs_agg_1w = nanmean(TW_oss_1w);
mean_mod_agg_1w = nanmean(TW_mod_1w);
covar_mod_agg_1w = nancov(TW_oss_1w,TW_mod_1w);
covar_mod_agg_1w = covar_mod_agg_1w(1,2);

alpha_agg_1w = std_mod_agg_1w/std_obs_agg_1w;
beta_agg_1w = mean_mod_agg_1w/mean_obs_agg_1w;
gamma_agg_1w = covar_mod_agg_1w/(std_obs_agg_1w*std_mod_agg_1w);
KGE_agg_1w = 1- sqrt((std_mod_agg_1w/std_obs_agg_1w-1)^2+(mean_mod_agg_1w/mean_obs_agg_1w-1)^2+...
    (covar_mod_agg_1w/(std_mod_agg_1w*std_obs_agg_1w)-1)^2);


sopra_agg_1w = [NaN(length(TW_oss_1w),1)];
sotto_agg_1w = [NaN(length(TW_oss_1w),1)];
sopra_agg_1w = nansum((TW_oss_1w - TW_mod_1w).^2);
sotto_agg_1w = nansum((TW_oss_1w - mean_obs_agg_1w).^2);
NSE_agg_1w = 1-(sopra_agg_1w/sotto_agg_1w);

% %----> errors at monthly time scale but with calibration at daily time
% scale

L1m = floor(length(Tw_oss_calib)/30);
TW_oss_1m = [NaN(L1m,1)];
TW_mod_1m = [NaN(L1m,1)];
p=0;
for i = 1:L1m
    TW_oss_1m(i)=nanmean(Tw_oss_calib(p+1:p+30));
    TW_mod_1m(i)=nanmean(Tw_sim_d_calib(p+1:p+30));
    p=p+30;
end

RMS_agg_1m = (nanmean( (TW_oss_1m-TW_mod_1m).^2 ))^0.5;

MAE_agg_1m = nanmean(abs(TW_oss_1m-TW_mod_1m));

std_obs_agg_1m = nanstd(TW_oss_1m);
std_mod_agg_1m = nanstd(TW_mod_1m);
mean_obs_agg_1m = nanmean(TW_oss_1m);
mean_mod_agg_1m = nanmean(TW_mod_1m);
covar_mod_agg_1m = nancov(TW_oss_1m,TW_mod_1m);
covar_mod_agg_1m = covar_mod_agg_1m(1,2);

alpha_agg_1m = std_mod_agg_1m/std_obs_agg_1m;
beta_agg_1m = mean_mod_agg_1m/mean_obs_agg_1m;
gamma_agg_1m = covar_mod_agg_1m/(std_obs_agg_1m*std_mod_agg_1m);
KGE_agg_1m = 1- sqrt((std_mod_agg_1m/std_obs_agg_1m-1)^2+(mean_mod_agg_1m/mean_obs_agg_1m-1)^2+...
    (covar_mod_agg_1m/(std_mod_agg_1m*std_obs_agg_1m)-1)^2);


sopra_agg_1m = [NaN(length(TW_oss_1m),1)];
sotto_agg_1m = [NaN(length(TW_oss_1m),1)];
sopra_agg_1m = nansum((TW_oss_1m - TW_mod_1m).^2);
sotto_agg_1m = nansum((TW_oss_1m - mean_obs_agg_1m).^2);
NSE_agg_1m = 1-(sopra_agg_1m/sotto_agg_1m);

% %----> errors at seasonal time scale but with calibration at daily time
% scale

prim_12w = isnan(Tw_oss_calib);
prim_12w = find(prim_12w==0);
prim_12w = prim_12w(1);
L12w = floor(length(Tw_oss_calib(prim_12w:end))/84);
TW_oss_12w = [NaN(L12w,1)];
TW_mod_12w = [NaN(L12w,1)];
p=prim_12w;
for i = 1:L12w
    TW_oss_12w(i)=nanmean(Tw_oss_calib(p+1:p+84));
    TW_mod_12w(i)=nanmean(Tw_sim_d_calib(p+1:p+84));
    p=p+84;
end

RMS_agg_12w = (nanmean( (TW_oss_12w-TW_mod_12w).^2 ))^0.5;

MAE_agg_12w = nanmean(abs(TW_oss_12w-TW_mod_12w));

std_obs_agg_12w = nanstd(TW_oss_12w);
std_mod_agg_12w = nanstd(TW_mod_12w);
mean_obs_agg_12w = nanmean(TW_oss_12w);
mean_mod_agg_12w = nanmean(TW_mod_12w);
covar_mod_agg_12w = nancov(TW_oss_12w,TW_mod_12w);
covar_mod_agg_12w = covar_mod_agg_12w(1,2);

alpha_agg_12w = std_mod_agg_12w/std_obs_agg_12w;
beta_agg_12w = mean_mod_agg_12w/mean_obs_agg_12w;
gamma_agg_12w = covar_mod_agg_12w/(std_obs_agg_12w*std_mod_agg_12w);
KGE_agg_12w = 1- sqrt((std_mod_agg_12w/std_obs_agg_12w-1)^2+(mean_mod_agg_12w/mean_obs_agg_12w-1)^2+...
    (covar_mod_agg_12w/(std_mod_agg_12w*std_obs_agg_12w)-1)^2);


sopra_agg_12w = [NaN(length(TW_oss_12w),1)];
sotto_agg_12w = [NaN(length(TW_oss_12w),1)];
sopra_agg_12w = nansum((TW_oss_12w - TW_mod_12w).^2);
sotto_agg_12w = nansum((TW_oss_12w - mean_obs_agg_12w).^2);
NSE_agg_12w = 1-(sopra_agg_12w/sotto_agg_12w);

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 1w
%----------------------------------------------------------------------

WW = strcmp(time_resolution,'1w');
if WW == 1
    Tw_1w_oss_cal = dati_cal(:,7);
    Tw_1w_mod_cal = dati_cal(:,8);
    
    RMS_cal_1w = (nanmean( (Tw_1w_oss_cal-Tw_1w_mod_cal).^2 ))^0.5;
    
    MAE_cal_1w = nanmean(abs(Tw_1w_oss_cal-Tw_1w_mod_cal));
    
    std_obs_cal_1w = nanstd( Tw_1w_oss_cal);
    std_mod_cal_1w = nanstd(Tw_1w_mod_cal);
    mean_obs_cal_1w = nanmean( Tw_1w_oss_cal);
    mean_mod_cal_1w = nanmean(Tw_1w_mod_cal);
    covar_mod_cal_1w = nancov( Tw_1w_oss_cal,Tw_1w_mod_cal);
    covar_mod_cal_1w = covar_mod_cal_1w(1,2);
    
    alpha_cal_1w = std_mod_cal_1w/std_obs_cal_1w;
    beta_cal_1w = mean_mod_cal_1w/mean_obs_cal_1w;
    gamma_cal_1w = covar_mod_cal_1w/(std_obs_cal_1w*std_mod_cal_1w);
    KGE_cal_1w = 1- sqrt((std_mod_cal_1w/std_obs_cal_1w-1)^2+(mean_mod_cal_1w/mean_obs_cal_1w-1)^2+...
        (covar_mod_cal_1w/(std_mod_cal_1w*std_obs_cal_1w)-1)^2);
    
    
    sopra_cal_1w = [NaN(length(Tw_1w_oss_cal),1)];
    sotto_cal_1w = [NaN(length(Tw_1w_oss_cal),1)];
    sopra_cal_1w = nansum((Tw_1w_oss_cal - Tw_1w_mod_cal).^2);
    sotto_cal_1w = nansum((Tw_1w_oss_cal - mean_obs_cal_1w).^2);
    NSE_cal_1w = 1-(sopra_cal_1w/sotto_cal_1w);
    
end

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 1m
%----------------------------------------------------------------------

MM = strcmp(time_resolution,'1m');
if MM == 1
    Tw_1m_oss_cal = dati_cal(:,7);
    Tw_1m_mod_cal = dati_cal(:,8);
    
    RMS_cal_1m = (nanmean( (Tw_1m_oss_cal-Tw_1m_mod_cal).^2 ))^0.5;
    
    MAE_cal_1m = nanmean(abs(Tw_1m_oss_cal-Tw_1m_mod_cal));
    
    std_obs_cal_1m = nanstd( Tw_1m_oss_cal);
    std_mod_cal_1m = nanstd(Tw_1m_mod_cal);
    mean_obs_cal_1m = nanmean( Tw_1m_oss_cal);
    mean_mod_cal_1m = nanmean(Tw_1m_mod_cal);
    covar_mod_cal_1m = nancov( Tw_1m_oss_cal,Tw_1m_mod_cal);
    covar_mod_cal_1m = covar_mod_cal_1m(1,2);
    
    alpha_cal_1m = std_mod_cal_1m/std_obs_cal_1m;
    beta_cal_1m = mean_mod_cal_1m/mean_obs_cal_1m;
    gamma_cal_1m = covar_mod_cal_1m/(std_obs_cal_1m*std_mod_cal_1m);
    KGE_cal_1m = 1- sqrt((std_mod_cal_1m/std_obs_cal_1m-1)^2+(mean_mod_cal_1m/mean_obs_cal_1m-1)^2+...
        (covar_mod_cal_1m/(std_mod_cal_1m*std_obs_cal_1m)-1)^2);
    
    
    sopra_cal_1m = [NaN(length(Tw_1m_oss_cal),1)];
    sotto_cal_1m = [NaN(length(Tw_1m_oss_cal),1)];
    sopra_cal_1m = nansum((Tw_1m_oss_cal - Tw_1m_mod_cal).^2);
    sotto_cal_1m = nansum((Tw_1m_oss_cal - mean_obs_cal_1m).^2);
    NSE_cal_1m = 1-(sopra_cal_1m/sotto_cal_1m);
    
end

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 12w
%----------------------------------------------------------------------

ST = strcmp(time_resolution,'12w');
if ST == 1
    Tw_12w_oss_cal = dati_cal(:,7);
    Tw_12w_mod_cal = dati_cal(:,8);
    
    RMS_cal_12w = (nanmean( (Tw_12w_oss_cal-Tw_12w_mod_cal).^2 ))^0.5;
    
    MAE_cal_12w = nanmean(abs(Tw_12w_oss_cal-Tw_12w_mod_cal));
    
    std_obs_cal_12w = nanstd( Tw_12w_oss_cal);
    std_mod_cal_12w = nanstd(Tw_12w_mod_cal);
    mean_obs_cal_12w = nanmean( Tw_12w_oss_cal);
    mean_mod_cal_12w = nanmean(Tw_12w_mod_cal);
    covar_mod_cal_12w = nancov( Tw_12w_oss_cal,Tw_12w_mod_cal);
    covar_mod_cal_12w = covar_mod_cal_12w(1,2);
    
    alpha_cal_12w = std_mod_cal_12w/std_obs_cal_12w;
    beta_cal_12w = mean_mod_cal_12w/mean_obs_cal_12w;
    gamma_cal_12w = covar_mod_cal_12w/(std_obs_cal_12w*std_mod_cal_12w);
    KGE_cal_12w = 1- sqrt((std_mod_cal_12w/std_obs_cal_12w-1)^2+(mean_mod_cal_12w/mean_obs_cal_12w-1)^2+...
        (covar_mod_cal_12w/(std_mod_cal_12w*std_obs_cal_12w)-1)^2);
    
    
    sopra_cal_12w = [NaN(length(Tw_12w_oss_cal),1)];
    sotto_cal_12w = [NaN(length(Tw_12w_oss_cal),1)];
    sopra_cal_12w = nansum((Tw_12w_oss_cal - Tw_12w_mod_cal).^2);
    sotto_cal_12w = nansum((Tw_12w_oss_cal - mean_obs_cal_12w).^2);
    NSE_cal_12w = 1-(sopra_cal_12w/sotto_cal_12w);
    
end

%----------------------------------------------------------------------
%---------------> ERRORS in validation period for daily data
%----------------------------------------------------------------------

RMS_val = (nanmean( (Tw_oss_valid-Tw_sim_d_valid).^2 ))^0.5;

MAE_val = nanmean(abs(Tw_oss_valid-Tw_sim_d_valid));

std_obs_val = nanstd(Tw_oss_valid);
std_mod_val = nanstd(Tw_sim_d_valid);
mean_obs_val = nanmean(Tw_oss_valid);
mean_mod_val = nanmean(Tw_sim_d_valid);
covar_mod_val = nancov(Tw_oss_valid,Tw_sim_d_valid);
covar_mod_val = covar_mod_val(1,2);

alpha_val = std_mod_val/std_obs_val;
beta_val = mean_mod_val/mean_obs_val;
gamma_val = covar_mod_val/(std_obs_val*std_mod_val);
KGE_val = 1- sqrt((std_mod_val/std_obs_val-1)^2+(mean_mod_val/mean_obs_val-1)^2+...
    (covar_mod_val/(std_mod_val*std_obs_val)-1)^2);


sopra_val = [NaN(length(Tw_oss_valid),1)];
sotto_val = [NaN(length(Tw_oss_valid),1)];
sopra_val = nansum((Tw_oss_valid - Tw_sim_d_valid).^2);
sotto_val = nansum((Tw_oss_valid - mean_obs_val).^2);
NSE_val = 1-(sopra_val/sotto_val);


% %----> errors at weekly time scale but with validation at daily time
% scale

L1w_val = floor(length(Tw_oss_valid)/7);
TW_oss_1w_val = [NaN(L1w_val,1)];
TW_mod_1w_val = [NaN(L1w_val,1)];
p=0;
for i = 1:L1w_val
    TW_oss_1w_val(i)=nanmean(Tw_oss_valid(p+1:p+7));
    TW_mod_1w_val(i)=nanmean(Tw_sim_d_valid(p+1:p+7));
    p=p+7;
end

RMS_agg_1w_val = (nanmean( (TW_oss_1w_val-TW_mod_1w_val).^2 ))^0.5;

MAE_agg_1w_val = nanmean(abs(TW_oss_1w_val-TW_mod_1w_val));

std_obs_agg_1w_val = nanstd(TW_oss_1w_val);
std_mod_agg_1w_val = nanstd(TW_mod_1w_val);
mean_obs_agg_1w_val = nanmean(TW_oss_1w_val);
mean_mod_agg_1w_val = nanmean(TW_mod_1w_val);
covar_mod_agg_1w_val = nancov(TW_oss_1w_val,TW_mod_1w_val);
covar_mod_agg_1w_val = covar_mod_agg_1w_val(1,2);

alpha_agg_1w_val = std_mod_agg_1w_val/std_obs_agg_1w_val;
beta_agg_1w_val = mean_mod_agg_1w_val/mean_obs_agg_1w_val;
gamma_agg_1w_val = covar_mod_agg_1w_val/(std_obs_agg_1w_val*std_mod_agg_1w_val);
KGE_agg_1w_val = 1- sqrt((std_mod_agg_1w_val/std_obs_agg_1w_val-1)^2+(mean_mod_agg_1w_val/mean_obs_agg_1w_val-1)^2+...
    (covar_mod_agg_1w_val/(std_mod_agg_1w_val*std_obs_agg_1w_val)-1)^2);


sopra_agg_1w_val = [NaN(length(TW_oss_1w_val),1)];
sotto_agg_1w_val = [NaN(length(TW_oss_1w_val),1)];
sopra_agg_1w_val = nansum((TW_oss_1w_val - TW_mod_1w_val).^2);
sotto_agg_1w_val = nansum((TW_oss_1w_val - mean_obs_agg_1w_val).^2);
NSE_agg_1w_val = 1-(sopra_agg_1w_val/sotto_agg_1w_val);

% %----> errors at monthly time scale but with validation at daily time
% scale

L1m_val = floor(length(Tw_oss_valid)/30);
TW_oss_1m_val = [NaN(L1m_val,1)];
TW_mod_1m_val = [NaN(L1m_val,1)];
p=0;
for i = 1:L1m_val
    TW_oss_1m_val(i)=nanmean(Tw_oss_valid(p+1:p+30));
    TW_mod_1m_val(i)=nanmean(Tw_sim_d_valid(p+1:p+30));
    p=p+30;
end

RMS_agg_1m_val = (nanmean( (TW_oss_1m_val-TW_mod_1m_val).^2 ))^0.5;

MAE_agg_1m_val = nanmean(abs(TW_oss_1m_val-TW_mod_1m_val));

std_obs_agg_1m_val = nanstd(TW_oss_1m_val);
std_mod_agg_1m_val = nanstd(TW_mod_1m_val);
mean_obs_agg_1m_val = nanmean(TW_oss_1m_val);
mean_mod_agg_1m_val = nanmean(TW_mod_1m_val);
covar_mod_agg_1m_val = nancov(TW_oss_1m_val,TW_mod_1m_val);
covar_mod_agg_1m_val = covar_mod_agg_1m_val(1,2);

alpha_agg_1m_val = std_mod_agg_1m_val/std_obs_agg_1m_val;
beta_agg_1m_val = mean_mod_agg_1m_val/mean_obs_agg_1m_val;
gamma_agg_1m_val = covar_mod_agg_1m_val/(std_obs_agg_1m_val*std_mod_agg_1m_val);
KGE_agg_1m_val = 1- sqrt((std_mod_agg_1m_val/std_obs_agg_1m_val-1)^2+(mean_mod_agg_1m_val/mean_obs_agg_1m_val-1)^2+...
    (covar_mod_agg_1m_val/(std_mod_agg_1m_val*std_obs_agg_1m_val)-1)^2);


sopra_agg_1m_val = [NaN(length(TW_oss_1m_val),1)];
sotto_agg_1m_val = [NaN(length(TW_oss_1m_val),1)];
sopra_agg_1m_val = nansum((TW_oss_1m_val - TW_mod_1m_val).^2);
sotto_agg_1m_val = nansum((TW_oss_1m_val - mean_obs_agg_1m_val).^2);
NSE_agg_1m_val = 1-(sopra_agg_1m_val/sotto_agg_1m_val);

% %----> errors at seasonal time scale but with calibration at daily time
% scale

L12w_val = floor(length(Tw_oss_valid)/84);
TW_oss_12w_val = [NaN(L12w_val,1)];
TW_mod_12w_val = [NaN(L12w_val,1)];
p=0;
for i = 1:L12w_val
    TW_oss_12w_val(i)=nanmean(Tw_oss_valid(p+1:p+84));
    TW_mod_12w_val(i)=nanmean(Tw_sim_d_valid(p+1:p+84));
    p=p+84;
end

RMS_agg_12w_val = (nanmean( (TW_oss_12w_val-TW_mod_12w_val).^2 ))^0.5;

MAE_agg_12w_val = nanmean(abs(TW_oss_12w_val-TW_mod_12w_val));

std_obs_agg_12w_val = nanstd(TW_oss_12w_val);
std_mod_agg_12w_val = nanstd(TW_mod_12w_val);
mean_obs_agg_12w_val = nanmean(TW_oss_12w_val);
mean_mod_agg_12w_val = nanmean(TW_mod_12w_val);
covar_mod_agg_12w_val = nancov(TW_oss_12w_val,TW_mod_12w_val);
covar_mod_agg_12w_val = covar_mod_agg_12w_val(1,2);

alpha_agg_12w_val = std_mod_agg_12w_val/std_obs_agg_12w_val;
beta_agg_12w_val = mean_mod_agg_12w_val/mean_obs_agg_12w_val;
gamma_agg_12w_val = covar_mod_agg_12w_val/(std_obs_agg_12w_val*std_mod_agg_12w_val);
KGE_agg_12w_val = 1- sqrt((std_mod_agg_12w_val/std_obs_agg_12w_val-1)^2+(mean_mod_agg_12w_val/mean_obs_agg_12w_val-1)^2+...
    (covar_mod_agg_12w_val/(std_mod_agg_12w_val*std_obs_agg_12w_val)-1)^2);


sopra_agg_12w_val = [NaN(length(TW_oss_12w_val),1)];
sotto_agg_12w_val = [NaN(length(TW_oss_12w_val),1)];
sopra_agg_12w_val = nansum((TW_oss_12w_val - TW_mod_12w_val).^2);
sotto_agg_12w_val = nansum((TW_oss_12w_val - mean_obs_agg_12w_val).^2);
NSE_agg_12w_val = 1-(sopra_agg_12w_val/sotto_agg_12w_val);

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 1w
%----------------------------------------------------------------------

WW = strcmp(time_resolution,'1w');
if WW == 1
    Tw_1w_oss_val = dati_val(:,7);
    Tw_1w_mod_val = dati_val(:,8);
    
    RMS_val_1w = (nanmean( (Tw_1w_oss_val-Tw_1w_mod_val).^2 ))^0.5;
    
    MAE_val_1w = nanmean(abs(Tw_1w_oss_val-Tw_1w_mod_val));
    
    std_obs_val_1w = nanstd( Tw_1w_oss_val);
    std_mod_val_1w = nanstd(Tw_1w_mod_val);
    mean_obs_val_1w = nanmean( Tw_1w_oss_val);
    mean_mod_val_1w = nanmean(Tw_1w_mod_val);
    covar_mod_val_1w = nancov( Tw_1w_oss_val,Tw_1w_mod_val);
    covar_mod_val_1w = covar_mod_val_1w(1,2);
    
    
    alpha_val_1w = std_mod_val_1w/std_obs_val_1w;
    beta_val_1w = mean_mod_val_1w/mean_obs_val_1w;
    gamma_val_1w = covar_mod_val_1w/(std_obs_val_1w*std_mod_val_1w);
    KGE_val_1w = 1- sqrt((std_mod_val_1w/std_obs_val_1w-1)^2+(mean_mod_val_1w/mean_obs_val_1w-1)^2+...
        (covar_mod_val_1w/(std_mod_val_1w*std_obs_val_1w)-1)^2);
    
    
    sopra_val_1w = [NaN(length(Tw_1w_oss_val),1)];
    sotto_val_1w = [NaN(length(Tw_1w_oss_val),1)];
    sopra_val_1w = nansum((Tw_1w_oss_val - Tw_1w_mod_val).^2);
    sotto_val_1w = nansum((Tw_1w_oss_val - mean_obs_val_1w).^2);
    NSE_val_1w = 1-(sopra_val_1w/sotto_val_1w);
    
end

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 1m
%----------------------------------------------------------------------

MM = strcmp(time_resolution,'1m');
if MM == 1
    Tw_1m_oss_val = dati_val(:,7);
    Tw_1m_mod_val = dati_val(:,8);
    
    RMS_val_1m = (nanmean( (Tw_1m_oss_val-Tw_1m_mod_val).^2 ))^0.5;
    
    MAE_val_1m = nanmean(abs(Tw_1m_oss_val-Tw_1m_mod_val));
    
    std_obs_val_1m = nanstd( Tw_1m_oss_val);
    std_mod_val_1m = nanstd(Tw_1m_mod_val);
    mean_obs_val_1m = nanmean( Tw_1m_oss_val);
    mean_mod_val_1m = nanmean(Tw_1m_mod_val);
    covar_mod_val_1m = nancov( Tw_1m_oss_val,Tw_1m_mod_val);
    covar_mod_val_1m = covar_mod_val_1m(1,2);
    
    alpha_val_1m = std_mod_val_1m/std_obs_val_1m;
    beta_val_1m = mean_mod_val_1m/mean_obs_val_1m;
    gamma_val_1m = covar_mod_val_1m/(std_obs_val_1m*std_mod_val_1m);
    KGE_val_1m = 1- sqrt((std_mod_val_1m/std_obs_val_1m-1)^2+(mean_mod_val_1m/mean_obs_val_1m-1)^2+...
        (covar_mod_val_1m/(std_mod_val_1m*std_obs_val_1m)-1)^2);
    
    
    sopra_val_1m = [NaN(length(Tw_1m_oss_val),1)];
    sotto_val_1m = [NaN(length(Tw_1m_oss_val),1)];
    sopra_val_1m = nansum((Tw_1m_oss_val - Tw_1m_mod_val).^2);
    sotto_val_1m = nansum((Tw_1m_oss_val - mean_obs_val_1m).^2);
    NSE_val_1m = 1-(sopra_val_1m/sotto_val_1m);
    
end

%----------------------------------------------------------------------
%---------------> ERRORS in calibration period for time_resolution = 12w
%----------------------------------------------------------------------

ST = strcmp(time_resolution,'12w');
if ST == 1
    Tw_12w_oss_val = dati_val(:,7);
    Tw_12w_mod_val = dati_val(:,8);
    
    RMS_val_12w = (nanmean( (Tw_12w_oss_val-Tw_12w_mod_val).^2 ))^0.5;
    
    MAE_val_12w = nanmean(abs(Tw_12w_oss_val-Tw_12w_mod_val));
    
    std_obs_val_12w = nanstd( Tw_12w_oss_val);
    std_mod_val_12w = nanstd(Tw_12w_mod_val);
    mean_obs_val_12w = nanmean( Tw_12w_oss_val);
    mean_mod_val_12w = nanmean(Tw_12w_mod_val);
    covar_mod_val_12w = nancov( Tw_12w_oss_val,Tw_12w_mod_val);
    covar_mod_val_12w = covar_mod_val_12w(1,2);
    
    alpha_val_12w = std_mod_val_12w/std_obs_val_12w;
    beta_val_12w = mean_mod_val_12w/mean_obs_val_12w;
    gamma_val_12w = covar_mod_val_12w/(std_obs_val_12w*std_mod_val_12w);
    KGE_val_12w = 1- sqrt((std_mod_val_12w/std_obs_val_12w-1)^2+(mean_mod_val_12w/mean_obs_val_12w-1)^2+...
        (covar_mod_val_12w/(std_mod_val_12w*std_obs_val_12w)-1)^2);
    
    
    sopra_val_12w = [NaN(length(Tw_12w_oss_val),1)];
    sotto_val_12w = [NaN(length(Tw_12w_oss_val),1)];
    sopra_val_12w = nansum((Tw_12w_oss_val - Tw_12w_mod_val).^2);
    sotto_val_12w = nansum((Tw_12w_oss_val - mean_obs_val_12w).^2);
    NSE_val_12w = 1-(sopra_val_12w/sotto_val_12w);
    
end
%-----------------------> saving errors
k=1;
ERRORI_cal(1,k)=RMS_cal;
ERRORI_cal(2,k)=MAE_cal;
ERRORI_cal(3,k)=KGE_cal;
ERRORI_cal(4,k)=NSE_cal;
ERRORI_cal(5,k)=alpha_cal;
ERRORI_cal(6,k)=beta_cal;
ERRORI_cal(7,k)=gamma_cal;

ERRORI_cal(8,k)=RMS_agg_1w;
ERRORI_cal(9,k)=MAE_agg_1w;
ERRORI_cal(10,k)=KGE_agg_1w;
ERRORI_cal(11,k)=NSE_agg_1w;
ERRORI_cal(12,k)=alpha_agg_1w;
ERRORI_cal(13,k)=beta_agg_1w;
ERRORI_cal(14,k)=gamma_agg_1w;

ERRORI_cal(15,k)=RMS_agg_1m;
ERRORI_cal(16,k)=MAE_agg_1m;
ERRORI_cal(17,k)=KGE_agg_1m;
ERRORI_cal(18,k)=NSE_agg_1m;
ERRORI_cal(19,k)=alpha_agg_1m;
ERRORI_cal(20,k)=beta_agg_1m;
ERRORI_cal(21,k)=gamma_agg_1m;

ERRORI_cal(22,k)=RMS_agg_12w;
ERRORI_cal(23,k)=MAE_agg_12w;
ERRORI_cal(24,k)=KGE_agg_12w;
ERRORI_cal(25,k)=NSE_agg_12w;
ERRORI_cal(26,k)=alpha_agg_12w;
ERRORI_cal(27,k)=beta_agg_12w;
ERRORI_cal(28,k)=gamma_agg_12w;

ERRORI_val(1,k)=RMS_val;
ERRORI_val(2,k)=MAE_val;
ERRORI_val(3,k)=KGE_val;
ERRORI_val(4,k)=NSE_val;
ERRORI_val(5,k)=alpha_val;
ERRORI_val(6,k)=beta_val;
ERRORI_val(7,k)=gamma_val;

ERRORI_val(8,k)=RMS_agg_1w_val;
ERRORI_val(9,k)=MAE_agg_1w_val;
ERRORI_val(10,k)=KGE_agg_1w_val;
ERRORI_val(11,k)=NSE_agg_1w_val;
ERRORI_val(12,k)=alpha_agg_1w_val;
ERRORI_val(13,k)=beta_agg_1w_val;
ERRORI_val(14,k)=gamma_agg_1w_val;

ERRORI_val(15,k)=RMS_agg_1m_val;
ERRORI_val(16,k)=MAE_agg_1m_val;
ERRORI_val(17,k)=KGE_agg_1m_val;
ERRORI_val(18,k)=NSE_agg_1m_val;
ERRORI_val(19,k)=alpha_agg_1m_val;
ERRORI_val(20,k)=beta_agg_1m_val;
ERRORI_val(21,k)=gamma_agg_1m_val;

ERRORI_val(22,k)=RMS_agg_12w_val;
ERRORI_val(23,k)=MAE_agg_12w_val;
ERRORI_val(24,k)=KGE_agg_12w_val;
ERRORI_val(25,k)=NSE_agg_12w_val;
ERRORI_val(26,k)=alpha_agg_12w_val;
ERRORI_val(27,k)=beta_agg_12w_val;
ERRORI_val(28,k)=gamma_agg_12w_val;

if WW==1
    ERRORI_cal(8,k)=RMS_cal_1w;
    ERRORI_cal(9,k)=MAE_cal_1w;
    ERRORI_cal(10,k)=KGE_cal_1w;
    ERRORI_cal(11,k)=NSE_cal_1w;
    ERRORI_cal(12,k)=alpha_cal_1w;
    ERRORI_cal(13,k)=beta_cal_1w;
    ERRORI_cal(14,k)=gamma_cal_1w;
    
    ERRORI_val(8,k)=RMS_val_1w;
    ERRORI_val(9,k)=MAE_val_1w;
    ERRORI_val(10,k)=KGE_val_1w;
    ERRORI_val(11,k)=NSE_val_1w;
    ERRORI_val(12,k)=alpha_val_1w;
    ERRORI_val(13,k)=beta_val_1w;
    ERRORI_val(14,k)=gamma_val_1w;
end

if MM==1
    ERRORI_cal(15,k)=RMS_cal_1m;
    ERRORI_cal(16,k)=MAE_cal_1m;
    ERRORI_cal(17,k)=KGE_cal_1m;
    ERRORI_cal(18,k)=NSE_cal_1m;
    ERRORI_cal(19,k)=alpha_cal_1m;
    ERRORI_cal(20,k)=beta_cal_1m;
    ERRORI_cal(21,k)=gamma_cal_1m;
    
    ERRORI_val(15,k)=RMS_val_1m;
    ERRORI_val(16,k)=MAE_val_1m;
    ERRORI_val(17,k)=KGE_val_1m;
    ERRORI_val(18,k)=NSE_val_1m;
    ERRORI_val(19,k)=alpha_val_1m;
    ERRORI_val(20,k)=beta_val_1m;
    ERRORI_val(21,k)=gamma_val_1m;
end

if ST==1
    ERRORI_cal(22,k)=RMS_cal_12w;
    ERRORI_cal(23,k)=MAE_cal_12w;
    ERRORI_cal(24,k)=KGE_cal_12w;
    ERRORI_cal(25,k)=NSE_cal_12w;
    ERRORI_cal(26,k)=alpha_cal_12w;
    ERRORI_cal(27,k)=beta_cal_12w;
    ERRORI_cal(28,k)=gamma_cal_12w;
    
    ERRORI_val(22,k)=RMS_val_12w;
    ERRORI_val(23,k)=MAE_val_12w;
    ERRORI_val(24,k)=KGE_val_12w;
    ERRORI_val(25,k)=NSE_val_12w;
    ERRORI_val(26,k)=alpha_val_12w;
    ERRORI_val(27,k)=beta_val_12w;
    ERRORI_val(28,k)=gamma_val_12w;
end

%----------------------------------------------------------------------
%------------------------> PLOTS
%----------------------------------------------------------------------

span=30;    % amplitude of the windoy which with the data are smoothed in the plots



if strcmp(time_resolution,'1d')                                      % only for time_resolution = daily
    figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 10]);
    plot(Ta_calib,'c')
    hold on
    plot(Tw_oss_calib,'b','linewidth',1)
    hold on
    plot(Tw_sim_d_calib,'r','linewidth',1)
    axis([0 length(Ta_calib) min(Ta_calib)-5 max(Ta_calib)+1]);
    name2 = {['Daily data in calibration period ' stazione_a '-' stazione_w ' ' time_resolution ...
        ' npar=' num2str(numero_parametri)]; [' KGE=' num2str(KGE_cal) ' RMS=' num2str(RMS_cal)...
        ' MAE=' num2str(MAE_cal) ' NSE=' num2str(NSE_cal)]};
    title(name2,'FontSize',15);
    xlabel('Time [d]','FontSize',15);
    ylabel('Temperature [°C]','FontSize',15);
    h_leg = legend('air temperature','water temperature',...
        'simulated water temperature');
    set(h_leg,'Location','SouthWest','FontSize',15)
    fnam1=['TA-TWoss-TWsim_d_calib_' stazione_a '_' stazione_w '_' time_resolution '.png'];
    print('-dpng','-r300', fullfile(fpat, fnam1));
    close all
    
    
    
    
    figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 10]);
    plot(nansmooth(Ta_calib,30),'c')
    hold on
    plot(nansmooth(Tw_oss_calib,30),'b','linewidth',1)
    hold on
    plot(nansmooth(Tw_sim_d_calib,30),'r','linewidth',1)
    axis([0 length(Ta_calib) min(Ta_calib)-5 max(Ta_calib)+1]);
    name10 = {['Smoothed daily data in calibration period ' stazione_a '-' stazione_w ' ' time_resolution ...
        ' npar=' num2str(numero_parametri)]; [' KGE=' num2str(KGE_cal) ' RMS=' num2str(RMS_cal)...
        ' MAE=' num2str(MAE_cal) ' NSE=' num2str(NSE_cal)]};
    title(name10,'FontSize',15);
    xlabel('Time [d]','FontSize',15);
    ylabel('Temperature [°C]','FontSize',15);
    h_leg = legend('air temperature','water temperature',...
        'simulated water temperature');
    set(h_leg,'Location','SouthWest','FontSize',15)
    fnam10=['TA-TWoss-TWsim_d_calib_smooth_' stazione_a '_' stazione_w '_' time_resolution '.png'];
    print('-dpng','-r300', fullfile(fpat, fnam10));
    close all
    
    
    
else
    
    figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 10]);
    plot(Ta_calib,'c')
    hold on
    plot(Tw_oss_calib,'b','linewidth',1)
    hold on
    plot(Tw_sim_d_calib,'r','linewidth',1)
    hold on
    plot(Tw_oss_aggr_calib,'ok')
    hold on
    plot(Tw_sim_aggr_calib,'+k')
    axis([0 length(Ta_calib) min(Ta_calib)-5 max(Ta_calib)+1]);
    name2 = {['Daily data in calibration period ' stazione_a '-' stazione_w ' ' time_resolution ...
        ' npar=' num2str(numero_parametri)]; [' KGE=' num2str(KGE_cal) ' RMS=' num2str(RMS_cal) ...
        ' MAE=' num2str(MAE_cal) ' NSE=' num2str(NSE_cal)]};
    title(name2,'FontSize',15);
    xlabel('Time [d]','FontSize',15);
    ylabel('Temperature [°C]','FontSize',15);
    h_leg = legend('air temperature','water temperature',...
        'simulated water temperature','water temperature','simulated water temperature');
    set(h_leg,'Location','SouthWest','FontSize',15)
    fnam1=['TA-TWoss-TWsim_d_calib_' stazione_a '_' stazione_w '_' time_resolution '.png'];
    print('-dpng','-r300', fullfile(fpat, fnam1));
    close all
end

figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 10]);
plot(Ta_valid,'c')
hold on
plot(Tw_oss_valid,'b','linewidth',1)
hold on
plot(Tw_sim_d_valid,'r','linewidth',1)
axis([0 length(Ta_valid) min(Ta_valid)-5 max(Ta_valid)+1]);
name3 = {['Daily data in validation period ' stazione_a '-' stazione_w ...
    ' ' time_resolution ' npar=' num2str(numero_parametri)]; [' KGE=' num2str(KGE_val)...
    ' RMS=' num2str(RMS_val) ' MAE=' num2str(MAE_val) ' NSE=' num2str(NSE_val)]};
title(name3,'FontSize',15);
xlabel('Time [d]','FontSize',15);
ylabel('Temperature [°C]','FontSize',15);
h_leg = legend('air temperature','water temperature',...
    'simulated water temperature');
set(h_leg,'Location','SouthWest','FontSize',15)
fnam2=['TA-TWoss-TWsim_d_valid_' stazione_a '_' stazione_w '_' time_resolution '.png'];
print('-dpng','-r300', fullfile(fpat, fnam2));
close all



figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 10]);
plot(nansmooth(Ta_valid,30),'c')
hold on
plot(nansmooth(Tw_oss_valid,30),'b','linewidth',1)
hold on
plot(nansmooth(Tw_sim_d_valid,30),'r','linewidth',1)
axis([0 length(Ta_valid) min(Ta_valid)-5 max(Ta_valid)+1]);
name30 = {['Smoothed daily data in validation period ' stazione_a '-' stazione_w ...
    ' ' time_resolution ' npar=' num2str(numero_parametri)]; [' KGE=' num2str(KGE_val)...
    ' RMS=' num2str(RMS_val) ' MAE=' num2str(MAE_val) ' NSE=' num2str(NSE_val)]};
title(name30,'FontSize',15);
xlabel('Time [d]','FontSize',15);
ylabel('Temperature [°C]','FontSize',15);
h_leg = legend('air temperature','water temperature',...
    'simulated water temperature');
set(h_leg,'Location','SouthWest','FontSize',15)
fnam20=['TA-TWoss-TWsim_d_valid_smooth_' stazione_a '_' stazione_w '_' time_resolution '.png'];
print('-dpng','-r300', fullfile(fpat, fnam20));
close all



figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 15 15]);
plot(smooth(TAam_calib,span), smooth(TWam_oss_calib,span), '.b')
grid on
hold on
plot(smooth(TAam_calib,span), smooth(TWam_sim_calib,span), '.r')
plot(bisettrice, bisettrice,'-k')
set(gca,'xlim',[Tminn_med,Tmaxx_med],'ylim',[Tminn_med,Tmaxx_med])
name4 = ['Hysteresis loop ' stazione_a '-' stazione_w ' ' time_resolution ...
    ' npar=' num2str(numero_parametri)];
title(name4,'FontSize',15);
xlabel('Air Temperature [°C]','FontSize',15);
ylabel('Water Temperature [°C]','FontSize',15);
h_leg =  legend('Ta-Tw.obs','Ta-Tw.mod');
set(h_leg,'Location','NorthWest','FontSize',15)
fnam3=['Isteresi_' stazione_a '_' stazione_w '_' time_resolution '.png'];
print('-dpng','-r300', fullfile(fpat, fnam3));
close all


figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 20]);
subplot(2,1,1)
plot(TAam_calib, 'c')
grid on
hold on
plot(TWam_oss_calib, 'b','linewidth',2)
plot(TWam_sim_calib, 'r','linewidth',2)
axis([0 366 min(TAam_calib)-10 max(TAam_calib)+2]);
name5 = {['Mean year in calibration period ' stazione_a '-' stazione_w ' ' time_resolution ' npar=' num2str(numero_parametri)];...
    [' KGE=' num2str(KGE_cal) ' RMS=' num2str(RMS_cal) ' MAE=' num2str(MAE_cal) ' NSE=' num2str(NSE_cal)]};
title(name5,'FontSize',15);
xlabel('Time [d]','FontSize',15);
ylabel('Temperature [°C]','FontSize',15);
h_leg = legend('air temperature am','water temperature am',...
    'simulated water temperature am');
set(h_leg,'Location','SouthWest','FontSize',15)
subplot(2,1,2)
plot(Qam_calib)
axis([0 366 0 max(Qam_calib)+1]);
grid on
xlabel('Time [d]','FontSize',15);
ylabel('Q [m^3/s]','FontSize',15);
fnam4=['TA-TWoss-TWsim_am_c_' stazione_a '_' stazione_w '_' time_resolution '.png'];
print('-dpng','-r300', fullfile(fpat, fnam4));
close all


figure('visible','off','PaperUnits','centimeters','PaperPosition',[1 1 20 20]);
subplot(2,1,1)
plot(TAam_valid, 'c')
grid on
hold on
plot(TWam_oss_valid, 'b','linewidth',2)
plot(TWam_sim_valid, 'r','linewidth',2)
axis([0 366 min(TAam_valid)-10 max(TAam_valid)+2]);
name58 = {['Mean year in validation period ' stazione_a '-' stazione_w ' ' time_resolution ' npar=' num2str(numero_parametri)];...
    [' KGE=' num2str(KGE_val) ' RMS=' num2str(RMS_val) ' MAE=' num2str(MAE_val) ' NSE=' num2str(NSE_val)]};
title(name58,'FontSize',15);
xlabel('Time [d]','FontSize',15);
ylabel('Temperature [°C]','FontSize',15);
h_leg = legend('air temperature am','water temperature am',...
    'simulated water temperature am');
set(h_leg,'Location','SouthWest','FontSize',15)
subplot(2,1,2)
plot(Qam_valid)
axis([0 366 min(Qam_valid)-1 max(Qam_valid)+1]);
grid on
xlabel('Time [d]','FontSize',15);
ylabel('Q [m^3/s]','FontSize',15);
fnam58=['TA-TWoss-TWsim_am_v_' stazione_a '_' stazione_w '_' time_resolution '.png'];
print('-dpng','-r300', fullfile(fpat, fnam58));
close all


fpat33=[nome_fiume '/errors/'];
file = fopen([fpat33 'errors' stazione_a '_' stazione_w '_npar' numero_parametri '_' time_resolution '.txt'],'w');
for k = 1:28
    fprintf(file,'%6.3f \t %6.3f \n',ERRORI_cal(k,1),ERRORI_val(k,1));
end
fclose(file);


