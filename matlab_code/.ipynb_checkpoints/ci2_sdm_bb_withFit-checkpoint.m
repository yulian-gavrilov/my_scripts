
omega_H_600 = 2.*pi*600.172821*1000000; % @600 MHz
omega_H_750 = 2.*pi*750.063525*1000000; % @750 MHz

measuredfreq = [0, 0, omega_H_600, omega_H_750, omega_H_600*2, omega_H_750*2];
measuredfreq_600old = [0, omega_H_600, omega_H_600*2];
measuredfreq_750old = [0, omega_H_750, omega_H_750*2];

gammah = 2.6752e8;
gamman = -2.71e7;

fieldH = 2*3.14159*[600.17,750.06]*1000000;
fieldN = -fieldH*gamman/gammah;

coeffH = 0.87; % 0.87
allfields = [0,0,fieldN, coeffH*fieldH];

measuredfreq = allfields;
measuredfreq_600 = [0,fieldN(1),coeffH*fieldH(1)];
measuredfreq_750 = [0,fieldN(2),coeffH*fieldH(2)];


sdfType = {'j0','jwx','jwh'};
field = {'600','750'};
System = {'WT','L49I','I57V','L49I_I57V'};
for k = 1:size(sdfType,2)
    for j = 1:size(field,2)
        for i = 1:size(System,2)

         eval (['load(''sdm_input/' sdfType{k} '_CI2_' System{i} '_' field{j} '.dat'')'])

        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 0*10^9:10^7:10*10^9;
%w = measuredfreq;
%%
Ser12_wt_s2_m3 = 0.77619642448060244; % m3-> m1 12S WT

Ser12_wt_s2_m7 = 0.6260092512516453;
Ser12_wt_s2f_m7 = 0.76341857106835209;
Ser12_wt_s2s_m7 = 0.82000788948006353;
Ser12_wt_ts_m7 = 6.8881218277587414e-09; % sec

Val13_wt_s2_m4 = 0.79873227127199065;
Val13_wt_ts_m4 = 2.9709565707888806e-12;

Val13_wt_s2_m7 = 0.79873227127199542;
Val13_wt_s2f_m7 = 1.0;
Val13_wt_s2s_m7 = 0.79873227127199542;
Val13_wt_ts_m7 = 2.9709565707849426e-12;

Lys11_wt_s2_m3 = 0.8165011164; % m3
%Lys11_wt_s2_m3 = 0.86121783; % m1


Lys11_wt_s2_m7 = 0.81650111;
Lys11_wt_s2f_m7 = 0.91652605;
Lys11_wt_s2s_m7 = 0.89086514;
Lys11_wt_ts_m7 = 4.148898e-24;
%%
Ser12_i57v_s2_m3 = 0.80879884353041176;

Ser12_i57v_s2_m4 = 0.79784536236718617;
Ser12_i57v_ts_m4 = 2.0056582580935269e-11;

Ser12_i57v_s2_m7 = 0.40651697503253692;
Ser12_i57v_s2f_m7 = 0.74519357560961186;
Ser12_i57v_s2s_m7 = 0.545518625412172;
Ser12_i57v_ts_m7 = 5.591318720144246e-09;
%%
Val13_l49i_s2_m4 = 0.8187908773394682;
Val13_l49i_ts_m4 = 1.8501922640187441e-11;

Val13_l49i_s2_m7 = 0.49703676110749068;
Val13_l49i_s2f_m7 = 0.79377232172838263;
Val13_l49i_s2s_m7 = 0.62617043641081938;
Val13_l49i_ts_m7 = 4.8513317287741962e-09;


Lys11_l49i_s2_m3 = 0.8284567; % m3
%Lys11_l49i_s2_m3 = 0.8668349; % m1

Lys11_l49i_s2_m7 = 0.8284567;
Lys11_l49i_s2f_m7 = 0.82863865;
Lys11_l49i_s2s_m7 = 0.999780463;
Lys11_l49i_ts_m7 = 3.39742784366e-27;
%%
tm = 3.62*10.^-9*1.2;
tm_wt = 3.44*10.^-9;
tm_l49i = 3.5906*10.^-9;
tm_i57v = 3.7898*10.^-9;

% tm_wt = tm;
% tm_l49i = tm;
% tm_i57v = tm;

for i = 1:size(w,2)
    
    Ser12_wt_m3(i) = jmf1(w(i),Ser12_wt_s2_m3,tm_wt);
    Ser12_wt_m7(i) = jmf3(w(i),Ser12_wt_s2_m7,Ser12_wt_s2f_m7,tm_wt,Ser12_wt_ts_m7);
    Ser12_i57v_m3(i) = jmf1(w(i),Ser12_i57v_s2_m3,tm_i57v);
    Ser12_i57v_m4(i) = jmf2(w(i),Ser12_i57v_s2_m4,tm_i57v,Ser12_i57v_ts_m4);
    Ser12_i57v_m7(i) = jmf3(w(i),Ser12_i57v_s2_m7,Ser12_i57v_s2f_m7,tm_i57v,Ser12_i57v_ts_m7);
    
    Lys11_wt_m3(i) = jmf1(w(i),Lys11_wt_s2_m3,tm_wt);
    Lys11_wt_m7(i) = jmf3(w(i),Lys11_wt_s2_m7,Lys11_wt_s2f_m7,tm_wt,Lys11_wt_ts_m7);

    Val13_wt_m4(i) = jmf2(w(i),Val13_wt_s2_m4,tm_wt,Val13_wt_ts_m4);
    Val13_wt_m7(i) = jmf3(w(i),Val13_wt_s2_m7,Val13_wt_s2f_m7,tm_wt,Val13_wt_ts_m7);
    Val13_l49i_m4(i) = jmf2(w(i),Val13_l49i_s2_m4,tm_l49i,Val13_l49i_ts_m4);
    Val13_l49i_m7(i) = jmf3(w(i),Val13_l49i_s2_m7,Val13_l49i_s2f_m7,tm_l49i,Val13_l49i_ts_m7);

    Lys11_l49i_m3(i) = jmf1(w(i),Lys11_l49i_s2_m3,tm_l49i);
    Lys11_l49i_m7(i) = jmf3(w(i),Lys11_l49i_s2_m7,Lys11_l49i_s2f_m7,tm_l49i,Lys11_l49i_ts_m7);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%ii=1;
%color1 = {'k','g','r','m'};
color1 = {'k','m'};

iii = [9,10,11];

for ii = 1:3
res = j0_CI2_WT_600(iii(ii),1);

if res == 12 
    System = {'WT','I57V'};
end

if res == 13 || res == 11
    System = {'WT','L49I'};
end

figure(ii)

hold on;
xlim ([-0.2*10^9, 5*10^9]);
A = max ([j0_CI2_WT_600(res,2),j0_CI2_WT_750(res,2),j0_CI2_L49I_600(res,2),j0_CI2_L49I_750(res,2),j0_CI2_I57V_600(res,2),j0_CI2_I57V_750(res,2),j0_CI2_L49I_I57V_600(res,2),j0_CI2_L49I_I57V_750(res,2)]);
B = A + 0.1*A;
C = -0.1*A;
ylim ([C, B]);

l = 0;
for i = 1:size(System,2)
    for k = 1:size(sdfType,2)
        for j = 1:size(field,2)
        
            eval (['p' num2str(l) ' = errorbar (measuredfreq_' field{j} '(' num2str(k) '),' sdfType{k} '_CI2_' System{i} '_' field{j} '(res,2),' sdfType{k} '_CI2_' System{i} '_' field{j} '(res,3),''MarkerSize'',5,''Marker'',''o'',''color'',''' color1{i} ''',''MarkerFaceColor'',''' color1{i} ''',''MarkerEdgeColor'',''none'',''LineStyle'',''none'');']);
            %eval (['p' num2str(l) ' = errorbar (measuredfreq_' field{j} '(' num2str(k) '),log10(' sdfType{k} '_CI2_' System{i} '_' field{j} '(res,2)), log10(' sdfType{k} '_CI2_' System{i} '_' field{j} '(res,3)),''MarkerSize'',5,''Marker'',''o'',''color'',''' color1{i} ''',''MarkerFaceColor'',''' color1{i} ''',''MarkerEdgeColor'',''none'',''LineStyle'',''none'');']);

            l = l+1;
        end
    end
end

x = w;
% if res == 12
%     pp1 = plot (x,log10(Ser12_wt_m3),'LineStyle','-','color','k');
%     pp2 = plot (x,log10(Ser12_wt_m7),'LineStyle','--','color','k');
%     
%     pp3 = plot (x,log10(Ser12_i57v_m3),'LineStyle','--','color','r');
%     pp4 = plot (x,log10(Ser12_i57v_m4),'LineStyle','--','color','r');
%     pp5 = plot (x,log10(Ser12_i57v_m7),'LineStyle','-','color','r');
%     leg=legend([p1, p6, pp1, pp2, pp3, pp4, pp5], 'WT', System{2}, 'WT - m3* (0.78)','WT - m7 (0.63 elim)','I57V - m3 (0.81)','I57V - m4 (0.80)','I57V - m7* (0.41)');
% 
% end
% 
% if res == 13
%     pp1 = plot (x,log10(Val13_wt_m4),'LineStyle','-','color','k');
%     pp2 = plot (x,log10(Val13_wt_m7),'LineStyle','--','color','k');
%     
%     pp3 = plot (x,log10(Val13_l49i_m4),'LineStyle','--','color','r');
%     pp4 = plot (x,log10(Val13_l49i_m7),'LineStyle','-','color','r');
%     leg=legend([p1, p6, pp1, pp2, pp3, pp4], 'WT', System{2}, 'WT - m4* (0.80)','WT - m7 (0.80)','L49I - m4 (0.82)','L49I - m7* (0.50)');
% 
% end

if res == 12
    pp1 = plot (x,Ser12_wt_m3,'LineStyle','-','color','k');
    pp2 = plot (x,Ser12_wt_m7,'LineStyle','--','color','k');
    
    pp3 = plot (x,Ser12_i57v_m3,'LineStyle','--','color','m');
    pp4 = plot (x,Ser12_i57v_m4,'LineStyle','--','color','m');
    pp5 = plot (x,Ser12_i57v_m7,'LineStyle','-','color','m');
    leg=legend([p1, p7, pp1, pp2, pp3, pp4, pp5], 'WT', System{2}, 'WT - m3* (0.78)','WT - m7 (0.63 elim)','I57V - m3 (0.81)','I57V - m4 (0.80)','I57V - m7* (0.41)');

end

if res == 13
    pp1 = plot (x,Val13_wt_m4,'LineStyle','-','color','k');
    pp2 = plot (x,Val13_wt_m7,'LineStyle','--','color','k');
    
    pp3 = plot (x,Val13_l49i_m4,'LineStyle','--','color','m');
    pp4 = plot (x,Val13_l49i_m7,'LineStyle','-','color','m');
    leg=legend([p1, p7, pp1, pp2, pp3, pp4], 'WT', System{2}, 'WT - m4* (0.80)','WT - m7 (0.80)','L49I - m4 (0.82)','L49I - m7* (0.50)');

end


if res == 11
    pp1 = plot (x,Lys11_wt_m3,'LineStyle','-','color','k');
    pp2 = plot (x,Lys11_wt_m7,'LineStyle','--','color','k');
    
    pp3 = plot (x,Lys11_l49i_m3,'LineStyle','-','color','m');
    pp4 = plot (x,Lys11_l49i_m7,'LineStyle','--','color','m');
    leg=legend([p1, p7, pp1, pp2, pp3, pp4], 'WT', System{2}, 'WT - m3* (0.82)','WT - m7 (0.82)','L49I - m3* (0.82)','L49I - m7 (0.82)');

end

%title ({'res#',num2str(res)});
title(sprintf('res #%s', num2str(res)),'FontName','nimbus roman no9 l','FontSize',25)                % Use ?sprintf?
set(gca,'FontName','nimbus roman no9 l','FontSize',20);
box on;
xlabel('\omega (s^-^1)','FontName','nimbus roman no9 l','FontSize',25);
ylabel('J (s)','FontName','nimbus roman no9 l','FontSize',25);

set(leg,'color','non','edgecolor','non','location','northeast')

eval(['print -dpng -r400 ''ci2_ ' num2str(res) '_SDM_fit.png''';]);


end


function jw1 = jmf1(w,s2,tm)
    jw1=(2.0/5.0)*(s2*tm/(1+w.^2*tm^2));
end

function jw2 = jmf2(w,s2,tm,ts)
    tsp2=tm*ts/(tm+ts);
    jw2=(2.0/5.0)*(s2*tm/(1+w.^2*tm^2)+(1-s2)*tsp2/(1+w.^2*tsp2.^2));
end

function jw3 = jmf3(w,s2,sf2,tm,ts)
    tsp3=tm*ts/(tm+ts);
    jw3=(2.0/5.0)*(s2*tm/(1+w.^2*tm^2)+(sf2-s2)*tsp3/(1+w.^2*tsp3.^2));
    %jw3=(2.0/5.0)*(s2*tm/(1+w.^2*tm^2)+sf2*tsp3/(1+w.^2*tsp3.^2));
end
