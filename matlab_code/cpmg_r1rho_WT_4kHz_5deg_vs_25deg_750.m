Data600 = importdata ('./cpmgallpeakser/ci2_wt_750MHZ_peaks_all.ser',' ',13);
Data750 = importdata ('./cpmgallpeakser/ci2_wt_5deg_750MHZ_peaks_all.ser',' ',13);

%vclist = load ('bb_relax_lists_cpmg_heroine.txt');
T = 0.05;

vclist1 = load ('bb_relax_lists_cpmg.txt');
vclist = vclist1./T;

load('./r1rho2and4kHz/R20_WT_R1r_750_4kHz');
load('./r1rho2and4kHz/R20_WT_5deg_R1r_750_4kHz');

numbRes = 59;

% vclist - in 50 ms
% x - in 1 s
% x = s*vclist/0.05 s
vclist_sec_600 = [vclist;1500];%./T;
vclist_sec_750 = [vclist;1500];%./T;


%mu_cpmg = 1./(2*vclist./1000); % in sec

% R2eff = -1/T*log(I/I_0)

for_title_600 = Data600.textdata(size (Data600.textdata,1)-numbRes+1:size(Data600.textdata,1),:);

for i = 1:size(Data600.data,1)

 resall = for_title_600 {i,7};    
 expression = '\d+';
 expression1 = '\w\d+';
 resNumb_600 = regexp(resall,expression,'match','once');    
 resName_600 = for_title_600{i,7}(1);
 resNameNumb_600 = regexp(resall,expression1,'match','once'); 
 
 resNameVec_600(i,1) = resName_600;
 resNumbVec_600(i,1) = str2num(resNumb_600);

%  resNameNumbVec_600{i,1} = resNumb_600;
  resNameNumbVec_600{i,1} = resNameNumb_600;


end

for_title_750 = Data750.textdata(size (Data750.textdata,1)-numbRes+1:size(Data750.textdata,1),:);

for i = 1:size(Data750.data,1)

 resall = for_title_750 {i,7};    
 expression = '\d+';
 expression1 = '\w\d+';
 resNumb_750 = regexp(resall,expression,'match','once');    
 resName_750 = for_title_750{i,7}(1); 
 resNameNumb_750 = regexp(resall,expression1,'match','once'); 
 
 resNameVec_750(i,1) = resName_750;
 resNumbVec_750(i,1) = str2num(resNumb_750);

%  resNameNumbVec_750{i,1} = resNumb_750;
  resNameNumbVec_750{i,1} = resNameNumb_750;

end

%%% TEMP !!!!!!!

%Data750.data(:,11:12)=NaN(size(Data750.data,1),2);

for i = 1:size(Data600.data,2)
    
     R2eff_600(:,i) = (-1/T)*log(Data600.data(:,i)./Data600.data(:,1));
     R2eff_750(:,i) = (-1/T)*log(Data750.data(:,i)./Data750.data(:,1));

     %R2eff_600(:,i) = (-1/T)*log(Data600.data(:,i)./Data600.data(:,1)); 
end

index1(:,1) = 1:1:size(Data600.data,1);

R2eff_600_forsort = horzcat (resNumbVec_600(:,1),index1,R2eff_600);
R2eff_750_forsort = horzcat (resNumbVec_750(:,1),index1,R2eff_750);

R2eff_600_sorted = sortrows(R2eff_600_forsort);
R2eff_750_sorted = sortrows(R2eff_750_forsort);

R2eff_600_sorted = horzcat(R2eff_600_sorted,R20_WT_R1r_750_4kHz(:,1));
R2eff_750_sorted = horzcat(R2eff_750_sorted,R20_WT_5deg_R1r_750_4kHz(:,1));

%a = resNameVec_600(R2eff_600(1));
j = 1;
tt = for_title_600{j,7}(1);


%figure('Renderer', 'painters', 'Position', [10 10 600 2000])
figure('Position', [10 10 6000 20000])

for plotId = 1:size(Data600.data,1)+1
    
    if plotId == size(Data600.data,1)+1
        
        h{plotId} = subplot(10,6,plotId);
        hLine{plotId} = plot(NaN,NaN,'color','k','LineStyle','none','MarkerFaceColor','k','Marker','o');
        hold on;
        plot(NaN,NaN,'color','r','LineStyle','none','MarkerFaceColor','r','Marker','o');
        axis off;
        h_legend=legend('@750MHz (WT at 25degC)','@750MHz  (WT at 5 degC)','Location','NorthEast');
        set(h_legend, 'Box', 'off','FontName','nimbus roman no9 l','FontSize',12);
        
    else
    
    res600 = resNameNumbVec_600{R2eff_600_sorted(plotId,2)};
    res750 = resNameNumbVec_750{R2eff_750_sorted(plotId,2)};
    res_600and750 = strcat(res600,' ',res750);
    
    x = vclist_sec_600;
    x1 = vclist_sec_750;    
    y = R2eff_600_sorted(plotId,3:end);
    y1 = R2eff_750_sorted(plotId,3:end);
    x (x==0)= NaN;
    
    h{plotId} = subplot(10,6,plotId);
    hLine{plotId} = plot(x, y,'MarkerSize',2,'Marker','o','color','k','MarkerFaceColor','k','MarkerEdgeColor','none','LineStyle','none');
    hold on;
    errorbar (1500, R20_WT_R1r_750_4kHz(plotId,1), R20_WT_R1r_750_4kHz(plotId,2),'MarkerSize',2,'Marker','o','color','k','MarkerFaceColor','k','MarkerEdgeColor','none','LineStyle','none')

    plot (x1, y1,'MarkerSize',2,'Marker','o','color','r','MarkerFaceColor','r','MarkerEdgeColor','none','LineStyle','none')
    
    errorbar (1500, R20_WT_5deg_R1r_750_4kHz(plotId,1), R20_WT_5deg_R1r_750_4kHz(plotId,2),'MarkerSize',2,'Marker','o','color','r','MarkerFaceColor','r','MarkerEdgeColor','none','LineStyle','none')

    %plot (x, y1,'MarkerSize',1,'Marker','o','color','k','MarkerFaceColor','k','MarkerEdgeColor','none','LineStyle','none')
    avgYY = mean([mean(y(:,20)),mean(y1(:,20))]);
    lineX = [1 4050];
    avgYYa = mean(y(:,20));
    avgYYb = mean(y1(:,20));
    avgYY2a = [avgYYa,avgYYa];
    avgYY2b = [avgYYb,avgYYb];
    %avgYY2 = [avgYY,avgYY];
    %line (lineX, avgYY2,'Color','k','LineStyle','--');
    
    line (lineX, avgYY2a,'Color','k','LineStyle','--');
    line (lineX, avgYY2b,'Color','r','LineStyle','--');

    hold on;
    box on;
    ylim ([0 35]);
    xlim ([1 1550]);

    
    xlabel('\nu_C_P_M_G (Hz)','FontName','nimbus roman no9 l','FontSize',5);
    ylabel('R_2_,_e_f_f','FontName','nimbus roman no9 l','FontSize',5);

    title(res600,'FontName','nimbus roman no9 l','FontSize',5);
    
    set(gca,'FontName','nimbus roman no9 l','FontSize',5);

    hold off;
    
    end
    
end

    print -dpng -r800 'ci2_wt_cpmg_r1rho_4kHz_5deg_vs_25deg_750.png';

