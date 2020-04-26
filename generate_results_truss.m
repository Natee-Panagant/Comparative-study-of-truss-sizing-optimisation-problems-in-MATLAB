clear all; close all; clc;
reset = 0; % reset =1 for re-calculate reference points, fronts and all metric results
recal = 1; % recal =1 for re-calculate all metric results
flag_conv = 1; % gather HV from all iteration for convergence plot


title_list={'10-bar';
            '25-bar';
            '37-bar';
            '60-bar';
            '72-bar';
            '120-bar';
            '200-bar';
            '942-bar'};
algo={  'MOALO';    %1
        'MODA';
        'MOGOA';
        'MOGWO';
        'MOMVO';
        'MOWCA';    %6
        'MSSA';
        'SHAMODE';
        'SHAMODE-WO';
        'NSGA-II';
        'RPBILDE';  %11
        'DEMO';
        'MOEA/D';
        'UPSEMOA'};
all_text=[];
load(['rst_' num2str(1,'%03.f') '_' num2str(1,'%03.f') '.mat']);
Niter=numel(rst(1).ppareto);
Nrun=numel(rst);
Ntest=numel(title_list);
Nalgo=numel(algo);
if reset==1
    % Gather References Points and Front
    RefPoint=cell(Ntest,1);
    RefPareto=cell(Ntest,1);
    RefFPareto=cell(Ntest,1);
    RefGPareto=cell(Ntest,1);
    A=cell(Ntest,1);
    for k=1:Ntest
        RefPoint{k}=zeros();
        RefPareto{k}=[];
        RefFPareto{k}=[];
        RefGPareto{k}=[];
        A=[];
        for i=1:Nalgo
            clc;
            ptext=['RefPoints gathering Progress = ' num2str(100*(((k-1)/Ntest)+(i/Nalgo)/Ntest),'%0.f') ' %%'];
            fprintf([all_text ptext]);
            fprintf('\n');
            for j=1:Nrun
                clear PPareto FPareto GPareto
                filename=['rst_' num2str(k,'%03.f') '_' num2str(i,'%03.f') '.mat'];
                load(filename);
                % Get reference points
                fea_ind=max(rst(j).gpareto{end},[],1)<=0;
                fi_max=max(rst(j).fpareto{end}(:,fea_ind),[],2);
                if numel(fi_max)==0
                    fi_max=zeros(size(rst(1).fpareto{1},1),1);
                end
                RefPoint{k}=max(RefPoint{k},fi_max);
                % Get reference front
                [RefPareto{k},RefFPareto{k},RefGPareto{k},A]=ref_front_sorting(rst(j).ppareto{end}(:,fea_ind),rst(j).fpareto{end}(:,fea_ind),rst(j).gpareto{end}(:,fea_ind),RefPareto{k},RefFPareto{k},RefGPareto{k},A,100);
            end
%             figure(1);hold on;
%             if i>1
%                 set(obj,'color','k');
%             end
%             obj=plot(RefFPareto{k}(1,:),RefFPareto{k}(2,:),'rx');
%             pause(0);
        end
    end
    ptext='RefPoints gathering Progress = *Complete';
    save RefPoint.mat RefPoint
    save RefFront.mat RefPareto RefFPareto RefGPareto
    all_text=[ptext '\n'];
    clc;
    fprintf(all_text);
end

if reset==1 || recal==1
    % Generate results
    result=struct();
    load RefPoint.mat
    load RefFront.mat

    % % Calculate All Metrics
    result.hv=zeros(Ntest,Nalgo,Nrun);
    result.hv_all_iter=zeros(Ntest,Nalgo,Nrun,Niter);
    result.gd=zeros(Ntest,Nalgo,Nrun);
    result.gd_all_iter=zeros(Ntest,Nalgo,Nrun,Niter);
    result.igd=zeros(Ntest,Nalgo,Nrun);
    result.igd_all_iter=zeros(Ntest,Nalgo,Nrun,Niter);
    result.ste=zeros(Ntest,Nalgo,Nrun);
    result.ste_all_iter=zeros(Ntest,Nalgo,Nrun,Niter);
    for k=1:Ntest
        for i=1:Nalgo
            clc;
            ptext=['Metrics calculation Progress = ' num2str(100*(((k-1)/Ntest)+(i/Nalgo)/Ntest),'%0.f') ' %%'];
            fprintf([all_text ptext]);
            %Load Results
            filename=['rst_' num2str(k,'%03.f') '_' num2str(i,'%03.f')];
            load(filename);
            for j=1:Nrun
                % Calculate All Metrics of all iteration for convergence plot
                if flag_conv==1
                    for iter=1:Niter
                        gpareto=rst(j).gpareto{iter};
                        fea_ind=max(gpareto,[],1)<=0; % extract only feasible solution
                        fpareto=rst(j).fpareto{iter}(:,fea_ind);
                        result.hv_all_iter(k,i,j,iter) = hypervolume(fpareto,RefPoint{k});
                        result.gd_all_iter(k,i,j,iter) = generational_distance(fpareto,RefFPareto{k});
                        result.igd_all_iter(k,i,j,iter) = inverted_generational_distance(fpareto,RefFPareto{k});
                        result.ste_all_iter(k,i,j,iter) = spacing_to_extent(fpareto);
                    end
                end
                % Calculate All Metrics of last iteration
                gpareto=rst(j).gpareto{end};
                fea_ind=max(rst(j).gpareto{end},[],1)<=0;
                result.hv(k,i,j) = hypervolume(rst(j).fpareto{end}(:,fea_ind),RefPoint{k});
                result.gd(k,i,j) = generational_distance(rst(j).fpareto{end}(:,fea_ind),RefFPareto{k});
                result.igd(k,i,j) = inverted_generational_distance(rst(j).fpareto{end}(:,fea_ind),RefFPareto{k});
                result.ste(k,i,j) = spacing_to_extent(rst(j).fpareto{end}(:,fea_ind));
            end
            [~,result.hv_ibest(k,i)]=max(reshape(result.hv_all_iter(k,i,:,Niter),1,[]));
        end
    end
    result.hv_mean=mean(result.hv,3);
    [result.hv_max,result.hv_imax]=max(result.hv,[],3);
    [result.hv_min,result.hv_imin]=min(result.hv,[],3);
    result.hv_std=std(result.hv,[],3);
    
    result.gd_mean=mean(result.gd,3);
    [result.gd_max,result.gd_imax]=max(result.gd,[],3);
    [result.gd_min,result.gd_imin]=min(result.gd,[],3);
    result.gd_std=std(result.gd,[],3);
    
    result.igd_mean=mean(result.igd,3);
    [result.igd_max,result.igd_imax]=max(result.igd,[],3);
    [result.igd_min,result.igd_imin]=min(result.igd,[],3);
    result.igd_std=std(result.igd,[],3);
    
    result.ste_mean=mean(result.ste,3);
    [result.ste_max,result.ste_imax]=max(result.ste,[],3);
    [result.ste_min,result.ste_imin]=min(result.ste,[],3);
    result.ste_std=std(result.ste,[],3);
    
    ptext='Metrics calculation Progress = *Complete';
    all_text=[all_text ptext '\n'];
    clc;
    fprintf(all_text);

    save('result.mat','result','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


load result.mat
for i=1:Ntest
    hv_i=reshape(result.hv(i,:,:),Nalgo,Nrun);
    [p,t,stats]=friedman(1./(hv_i'),1,'off');
    result.hv_rank(i,:)=stats.meanranks;
    
    gd_i=reshape(result.gd(i,:,:),Nalgo,Nrun);
    [p,t,stats]=friedman((gd_i'),1,'off');
    result.gd_rank(i,:)=stats.meanranks;
    
    igd_i=reshape(result.igd(i,:,:),Nalgo,Nrun);
    [p,t,stats]=friedman((igd_i'),1,'off');
    result.igd_rank(i,:)=stats.meanranks;
    
    ste_i=reshape(result.ste(i,:,:),Nalgo,Nrun);
    [p,t,stats]=friedman((ste_i'),1,'off');
    result.ste_rank(i,:)=stats.meanranks;
end

hv_table=[];
gd_table=[];
igd_table=[];
ste_table=[];

tab_hv=[];
tab_gd=[];
tab_igd=[];
tab_ste=[];
for i=1:Nalgo
    tab_hv=[tab_hv;result.hv_mean(:,i)'];
    tab_hv=[tab_hv;result.hv_max(:,i)'];
    tab_hv=[tab_hv;result.hv_min(:,i)'];
    tab_hv=[tab_hv;result.hv_std(:,i)'];
    tab_hv=[tab_hv;result.hv_rank(:,i)'];
    
    tab_gd=[tab_gd;result.gd_mean(:,i)'];
    tab_gd=[tab_gd;result.gd_max(:,i)'];
    tab_gd=[tab_gd;result.gd_min(:,i)'];
    tab_gd=[tab_gd;result.gd_std(:,i)'];
    tab_gd=[tab_gd;result.gd_rank(:,i)'];
    
    tab_igd=[tab_igd;result.igd_mean(:,i)'];
    tab_igd=[tab_igd;result.igd_max(:,i)'];
    tab_igd=[tab_igd;result.igd_min(:,i)'];
    tab_igd=[tab_igd;result.igd_std(:,i)'];
    tab_igd=[tab_igd;result.igd_rank(:,i)'];
    
    tab_ste=[tab_ste;result.ste_mean(:,i)'];
    tab_ste=[tab_ste;result.ste_max(:,i)'];
    tab_ste=[tab_ste;result.ste_min(:,i)'];
    tab_ste=[tab_ste;result.ste_std(:,i)'];
    tab_ste=[tab_ste;result.ste_rank(:,i)'];
end

hv_mean_best_algo=cell(Ntest,1);
hv_rank_best_algo=cell(Ntest,1);
gd_mean_best_algo=cell(Ntest,1);
gd_rank_best_algo=cell(Ntest,1);
igd_mean_best_algo=cell(Ntest,1);
igd_rank_best_algo=cell(Ntest,1);
ste_mean_best_algo=cell(Ntest,1);
ste_rank_best_algo=cell(Ntest,1);
result.ibest_mean_hv=zeros(Ntest,2);
result.ibest_mean_gd=zeros(Ntest,2);
result.ibest_mean_igd=zeros(Ntest,2);
result.ibest_mean_ste=zeros(Ntest,2);
result.ibest_rank_hv=zeros(Ntest,2);
result.ibest_rank_gd=zeros(Ntest,2);
result.ibest_rank_igd=zeros(Ntest,2);
result.ibest_rank_ste=zeros(Ntest,2);
for i=1:Ntest
    [~,result.ibest_mean_hv(i,1)]=max(result.hv_mean(i,:));
    [~,result.ibest_mean_hv(i,2)]=max(result.hv(i,result.ibest_mean_hv(i,1),:));
    hv_mean_best_algo{i,1}=algo{result.ibest_mean_hv(i,1)};
    [~,result.ibest_rank_hv(i,1)]=min(result.hv_rank(i,:));
    [~,result.ibest_rank_hv(i,2)]=min(result.hv(i,result.ibest_rank_hv(i,1),:));
    hv_rank_best_algo{i,1}=algo{result.ibest_rank_hv(i,1)};
    
    [~,result.ibest_mean_gd(i,1)]=min(result.gd_mean(i,:));
    [~,result.ibest_mean_gd(i,2)]=min(result.gd(i,result.ibest_mean_gd(i,1),:));
    gd_mean_best_algo{i,1}=algo{result.ibest_mean_gd(i,1)};
    [~,result.ibest_rank_gd(i,1)]=min(result.gd_rank(i,:));
    [~,result.ibest_rank_gd(i,2)]=min(result.gd(i,result.ibest_rank_gd(i,1),:));
    gd_rank_best_algo{i,1}=algo{result.ibest_rank_gd(i,1)};
    
    [~,result.ibest_mean_igd(i,1)]=min(result.igd_mean(i,:));
    [~,result.ibest_mean_igd(i,2)]=min(result.igd(i,result.ibest_mean_igd(i,1),:));
    igd_mean_best_algo{i,1}=algo{result.ibest_mean_igd(i,1)};
    [~,result.ibest_rank_igd(i,1)]=min(result.igd_rank(i,:));
    [~,result.ibest_rank_igd(i,2)]=min(result.igd(i,result.ibest_rank_igd(i,1),:));
    igd_rank_best_algo{i,1}=algo{result.ibest_rank_igd(i,1)};
    
    [~,result.ibest_mean_ste(i,1)]=min(result.ste_mean(i,:));
    [~,result.ibest_mean_ste(i,2)]=min(result.ste(i,result.ibest_mean_ste(i,1),:));
    ste_mean_best_algo{i,1}=algo{result.ibest_mean_ste(i,1)};
    [~,result.ibest_rank_ste(i,1)]=min(result.ste_rank(i,:));
    [~,result.ibest_rank_ste(i,2)]=min(result.ste(i,result.ibest_rank_ste(i,1),:));
    ste_rank_best_algo{i,1}=algo{result.ibest_rank_ste(i,1)};
end
result.overall_rank=((result.hv_rank+result.gd_rank+result.igd_rank+result.ste_rank)/4)';
[~,result.ibest_overall_rank]=min(result.overall_rank,[],1);
tab_overall_rank=result.overall_rank;
save('result.mat','result','tab_hv','tab_gd','tab_igd','tab_ste','tab_overall_rank');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tab_xx variable contents (xx = each metric) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size = [5(mean,max,min,std,Friedman_Rank)xNalgo] x Ntest
% For benchmark test -> Size = (5x14) x 21 = 70x21
tab_hv
tab_gd
tab_igd
tab_ste 
% algo1_test1_mean      algo1_test2_mean        . . .   algo1_test8_mean
% algo1_test1_max       algo1_test2_max         . . .   algo1_test8_max
% algo1_test1_min       algo1_test2_min         . . .   algo1_test8_min
% algo1_test1_std       algo1_test2_std         . . .   algo1_test8_std
% algo1_test1_Frank     algo1_test2_Frank       . . .   algo1_test8_Frank
%
% algo2_test1_mean      algo2_test2_mean        . . .   algo2_test8_mean
% algo2_test1_max       algo2_test2_max         . . .   algo2_test8_max
% algo2_test1_min       algo2_test2_min         . . .   algo2_test8_min
% algo2_test1_std       algo2_test2_std         . . .   algo2_test8_std
% algo2_test1_Frank     algo2_test2_Frank       . . .   algo2_test8_Frank
%       .                       .               .               .
%       .                       .                 .             .
%       .                       .                   .           .
% algo14_test1_mean     algo14_test2_mean       . . .   algo14_test8_mean
% algo14_test1_max      algo14_test2_max        . . .   algo14_test8_max
% algo14_test1_min      algo14_test2_min        . . .   algo14_test8_min
% algo14_test1_std      algo14_test2_std        . . .   algo14_test8_std
% algo14_test1_Frank    algo14_test2_Frank      . . .   algo14_test8_Frank

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tab_overall_rank variable contents (mean of Friedman rank of all metrics) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tab_overall_rank

% algo1_test1_mean_Frank	algo1_test2_mean_Frank      . . . 	algo1_test8_mean_Frank 
% algo2_test1_mean_Frank    algo2_test2_mean_Frank      . . .   algo2_test8_mean_Frank
% algo3_test1_mean_Frank    algo3_test2_mean_Frank      . . .   algo3_test8_mean_Frank 
%       .                       .                       .                   .
%       .                       .                         .                 .
%       .                       .                           .               .
% algo14_test1_mean_Frank    algo14_test2_mean_Frank    . . .   algo14_test8_mean_Frank 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc={'-xr'%
    '-xg'
    '-xb'
    '-xm'
    '-or'%
    '-og'
    '-ob'
    '-om'
    '-^r'%
    '-^g'
    '-^b'
    '-^m'
    '-sr'%
    '-sg'
    '-sb'
    '-sm'};
load result.mat
for m=1:4
    rs_conv=zeros(Nalgo,Niter);
    for k=1:Ntest
        figure(1);clf; hold on;
        xlo=xlabel('Iteration');
        switch m
            case 1
                rs=result.hv_all_iter;
                ylo=ylabel('Hypervolume');
                Cfilename=['hv_conv' num2str(k) '.tif'];
            case 2
                rs=result.gd_all_iter;
                ylo=ylabel('Generational Distance');
                Cfilename=['gd_conv' num2str(k) '.tif'];
            case 3
                rs=result.igd_all_iter;
                ylo=ylabel('Inverted Generational Distance');
                Cfilename=['igd_conv' num2str(k) '.tif'];
            case 4
                rs=result.ste_all_iter;
                ylo=ylabel('Spacing-to-Extent');
                Cfilename=['ste_conv' num2str(k) '.tif'];
        end
        set(xlo,'fontweight','bold');
        set(ylo,'fontweight','bold');
        for i=1:Nalgo
            rs_ki=reshape(mean(rs(k,i,:,:),3),1,Niter,1,1);
            rs_conv(i,:)=rs_ki;
            plot(1:Niter,rs_conv(i,:),mc{i});
        end
        lo=legend({'MOALO','MODA','MOGOA','MOGWO','MOMVO','MOWCA','MSSA','SHAMODE','SHAMODE-WO','NSGA-II',...
                'RPBILDE','DEMO','MOEA/D','UPSEMOA'}','Location','southeast','NumColumns',7,...
                'Orientation','horizontal','fontsize',10);
        set(gca,'unit','norm','position',[0.1 0.2 0.8 0.75],'fontsize',12); 
        set(gcf,'position',[100 100 900 600]);    
        set(lo,'position',[0.05,0.05,0.9,0.05]);
        title(title_list{k},'fontsize',15,'fontweight','bold');
        saveas(gcf,Cfilename);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best Pareto plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mc={'xr'%
    'xg'
    'xb'
    'xm'
    'or'%
    'og'
    'ob'
    'om'
    '^r'%
    '^g'
    '^b'
    '^m'
    'sr'%
    'sg'
    'sb'
    'sm'};
load RefPoint.mat
load RefFront.mat
for m=1:2
    rs_conv=zeros(Nalgo,Niter);
    for k=1:Ntest
        figure(1);clf; hold on;
        switch m
            case 1
                hvi(:,:)=result.hv(k,:,:);
                [~,ibest]=max(hvi,[],2);
                Bfilename=['best_HV_Pareto' num2str(k)];
                metric='HV';
            case 2
                igdi(:,:)=result.igd(k,:,:);
                [~,ibest]=max(igdi,[],2);
                Bfilename=['best_IGD_Pareto' num2str(k)];
                metric='IGD';
        end
        
        for i=1:Nalgo
            filename=['rst_' num2str(k,'%03.f') '_' num2str(i,'%03.f') '.mat'];
            load(filename);
            F=rst(ibest(i)).fpareto{end};
            [~,isort]=sort(F(1,:));
            F=F(:,isort);
            pobj=plot(F(1,:),F(2,:),mc{i});
            set(pobj,'markersize',4);
        end
        [~,sind]=sort(RefFPareto{k}(1,:));
        RefFPareto{k}=RefFPareto{k}(:,sind);
        robj=plot(RefFPareto{k}(1,:),RefFPareto{k}(2,:),'k-','linewidth',1);
        
        xlabel('Mass','fontweight','bold');
        ylabel('Compliance','fontweight','bold');
        title({['Best ' metric ' of ' title_list{k}]},'fontsize',15);
        set(gca,'fontsize',16)
        lo=legend({'MOALO','MODA','MOGOA','MOGWO','MOMVO','MOWCA','MSSA','SHAMODE','SHAMODE-WO','NSGA-II',...
                'RPBILDE','DEMO','MOEA/D','UPSEMOA','Reference Pareto'}','Location','southeast','NumColumns',2,...
                'Orientation','horizontal','fontsize',10);
        set(lo,'position',[0.45,0.58,0.4,0.3]);
        axis tight;
        print(Bfilename,'-dtiffn','-r300');
    end
end