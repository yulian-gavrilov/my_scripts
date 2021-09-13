tic;
% This script average different runnings of t48@94 (ubiquitinated ubc7)
% made with different force fields (mber99, amber99sb)
coord = importdata('coord_CUT6_md100ns1_75ns.txt');
coord1 = importdata('coord_CUT6_md100ns2_75ns.txt');%coord_CUT4_md100ns3.txt
coord2 = importdata('coord_CUT6_md100ns3_75ns.txt');%coord_CUT6_md100ns3.txt

% load('subWTavg');
% load('subCUT4avg');
% load('subCUT6avg');
% sub1 = subCUT4avg - subWTavg;
% sub2 = subCUT6avg - subWTavg;

distance = zeros(93, 93);
distance1 = zeros(93, 93);
distance2 = zeros(93, 93);



sum = zeros(93, 93); sum1 = zeros(93, 93); sum2 = zeros(93, 93);
%sum3D = zeros(93, 93, 10);
number=751; %input('what is the nummber of frames? ');
x = 0;
y = 0;
 for n=0:(number-1)
    for i = 1:93;
        x = 93*n + i;
            for j= 1:93;
                y = 93*n + j;
                distance(i,j) = ((coord(x,1)-coord(y,1))^2 + (coord(x,2)-coord(y,2))^2 + (coord(x,3)-coord(y,3))^2)^0.5;
                distance1(i,j) = ((coord1(x,1)-coord1(y,1))^2 + (coord1(x,2)-coord1(y,2))^2 + (coord1(x,3)-coord1(y,3))^2)^0.5;
                distance2(i,j) = ((coord2(x,1)-coord2(y,1))^2 + (coord2(x,2)-coord2(y,2))^2 + (coord2(x,3)-coord2(y,3))^2)^0.5;
            end
    end
    sum = sum + distance;
    sum1 = sum1 + distance1;
    sum2 = sum2 + distance2;
    sum3D{n+1} = distance;
    sum3D1{n+1} = distance1;
    sum3D2{n+1} = distance2;
 end

current_matrix= sum/(number);
current_matrix1= sum1/(number);
current_matrix2= sum2/(number);
average = (current_matrix + current_matrix1 + current_matrix2)./3;



  step = 10;


B = cell2mat(arrayfun(@(x)permute(x{:},[1 3 2]),sum3D,'UniformOutput',false));
B1 = cell2mat(arrayfun(@(x)permute(x{:},[1 3 2]),sum3D1,'UniformOutput',false));
B2 = cell2mat(arrayfun(@(x)permute(x{:},[1 3 2]),sum3D2,'UniformOutput',false));

%B11(:,:) = B1 (:,1,:);

s = std(B,0,2);s1 = std(B1,0,2);s2 = std(B2,0,2);

ss(:,:) = s(:,1,:);
ss1(:,:) = s1(:,1,:);
ss2(:,:) = s2(:,1,:);
%average_std = (ss + ss1 + ss2)./3;
average_std = sqrt(ss.^2 + ss1.^2 + ss2.^2);


    figure(1);
    h = pcolor(ss);
    set(h,'LineStyle','none');
    xl='Residue number';
    yl='Residue number';
    xlabel(xl,'FontSize', 14);
    ylabel(yl,'FontSize', 14);
    colorbar();%%caxis([0 50]);
    %caxis([0 6]);
    %set(gca,'FontSize',20);

axes_new0 = 1:93;

for i =1:93
if axes_new0(i) > 70
axes_new0(i) = axes_new0(i)+6;
end
end
 axes_new = axes_new0';

%axes1 = axes('Parent',figure1,'XTickLabel',axes_new,'FontSize',20);

set(gca,'XTick',1:step:93);
set(gca,'YTick',1:step:93);
set(gca,'XTickLabel',axes_new(1:step:end));
set(gca,'YTickLabel',axes_new(1:step:end));


    title ('Distance Matrice ACP CUT6. MD 1. Colorbar - distance in angstoms','FontSize', 15);
     %%print -dpng -r200 'DM_ACP_CUT6_MD1_std.png'

    figure(2);
    h = pcolor(ss1);
    set(h,'LineStyle','none');
    xl='Residue number';
    yl='Residue number';
    xlabel(xl,'FontSize', 14);
    ylabel(yl,'FontSize', 14);
    colorbar();%%caxis([0 50]);
    %%caxis([0 6]);
    %set(gca,'FontSize',20);

axes_new0 = 1:93;

for i =1:93
if axes_new0(i) > 70
axes_new0(i) = axes_new0(i)+6;
end
end
 axes_new = axes_new0';

%axes1 = axes('Parent',figure1,'XTickLabel',axes_new,'FontSize',20);

set(gca,'XTick',1:step:93);
set(gca,'YTick',1:step:93);
set(gca,'XTickLabel',axes_new(1:step:end));
set(gca,'YTickLabel',axes_new(1:step:end));


    title ('Distance Matrice ACP CUT6. MD 2. Colorbar - distance in angstoms','FontSize', 15);
     %%print -dpng -r200 'DM_ACP_CUT6_MD2_std.png'


    figure(3);
    h = pcolor(ss2);
    set(h,'LineStyle','none');
    xl='Residue number';
    yl='Residue number';
    xlabel(xl,'FontSize', 14);
    ylabel(yl,'FontSize', 14);
    colorbar();%%caxis([0 50]);
    %caxis([0 6]);
    %set(gca,'FontSize',20);

axes_new0 = 1:93;

for i =1:93
if axes_new0(i) > 70
axes_new0(i) = axes_new0(i)+6;
end
end
 axes_new = axes_new0';

%axes1 = axes('Parent',figure1,'XTickLabel',axes_new,'FontSize',20);

set(gca,'XTick',1:step:93);
set(gca,'YTick',1:step:93);
set(gca,'XTickLabel',axes_new(1:step:end));
set(gca,'YTickLabel',axes_new(1:step:end));


    title ('Distance Matrice ACP CUT6. MD 3. Colorbar - distance in angstoms','FontSize', 15);
     %%print -dpng -r200 'DM_ACP_CUT6_MD3_std.png'









%     figure(1);
%     h = pcolor(current_matrix);
%     set(h,'LineStyle','none');
%     xl='Residue number';
%     yl='Residue number';
%     xlabel(xl,'FontSize', 14);
%     ylabel(yl,'FontSize', 14);
%     colorbar();%%caxis([0 50]);
%     %caxis([0 50]);
%     %set(gca,'FontSize',20);
%
% axes_new0 = 1:93;
%
% for i =1:93
% if axes_new0(i) > 70
% axes_new0(i) = axes_new0(i)+6;
% end
% end
%  axes_new = axes_new0';
%
% %axes1 = axes('Parent',figure1,'XTickLabel',axes_new,'FontSize',20);
%
% set(gca,'XTick',1:step:93);
% set(gca,'YTick',1:step:93);
% set(gca,'XTickLabel',axes_new(1:step:end));
% set(gca,'YTickLabel',axes_new(1:step:end));
%
%
%     title ('Distance Matrice ACP WT. MD 3. Colorbar - distance in angstoms','FontSize', 15);
%     print -dpng -r200 'DM_ACP_WT_MD3.png'
%
%
%
%
%     figure(2);
%     h = pcolor(current_matrix1);
%     set(h,'LineStyle','none');
%     xl='Residue number';
%     yl='Residue number';
%     xlabel(xl,'FontSize', 14);
%     ylabel(yl,'FontSize', 14);
%     colorbar();%%caxis([0 50]);
%     %caxis([0 50]);
%     %set(gca,'FontSize',20);
%
% set(gca,'XTick',1:step:93);
% set(gca,'YTick',1:step:93);
% set(gca,'XTickLabel',axes_new(1:step:end));
% set(gca,'YTickLabel',axes_new(1:step:end));
%
%
%     title ('Distance Matrice ACP CUT4. MD 3. Colorbar - distance in angstoms','FontSize', 15);
%     print -dpng -r200 'DM_ACP_CUT4_MD3.png'
%
%     figure(3);
%     h = pcolor(current_matrix2);
%     set(h,'LineStyle','none');
%     xl='Residue number';
%     yl='Residue number';
%     xlabel(xl,'FontSize', 14);
%     ylabel(yl,'FontSize', 14);
%     colorbar();%%caxis([0 50]);
%     %caxis([0 50]);
%     %set(gca,'FontSize',20);
%
% set(gca,'XTick',1:step:93);
% set(gca,'YTick',1:step:93);
% set(gca,'XTickLabel',axes_new(1:step:end));
% set(gca,'YTickLabel',axes_new(1:step:end));
%
%
%     title ('Distance Matrice ACP CUT6. MD 3. Colorbar - distance in angstoms','FontSize', 15);
%     print -dpng -r200 'DM_ACP_CUT6_MD3.png'


    toc;
