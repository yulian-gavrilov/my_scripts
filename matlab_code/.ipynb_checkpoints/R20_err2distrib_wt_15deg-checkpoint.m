
%%%%%%%%

Data600 = importdata ('./cpmgallpeakser/ci2_wt_15deg_600MHZ_peaks_all.ser',' ',13); % unsorted res58, res30 etc
Data750 = importdata ('./cpmgallpeakser/ci2_wt_15deg_750MHZ_peaks_all.ser',' ',13); % unsorted

Data600_max = importdata('./serfiles_max_min/WT_15deg_600MHz_peaks20.ser');% unsorted
Data750_max = importdata('./serfiles_max_min/WT_15deg_750MHz_peaks20.ser');% unsorted

% load('./r1rho2and4kHz/R20_WT_R1r_600_2kHz');R20_600_2kHz = R20_WT_R1r_600_2kHz;% sorted! res 2, res3 etc
% load('./r1rho2and4kHz/R20_WT_R1r_750_2kHz');R20_750_2kHz = R20_WT_R1r_750_2kHz;% sorted!

load('./r1rho2and4kHz/R20_WT_15deg_R1r_600_4kHz');R20_600_4kHz = R20_WT_15deg_R1r_600_4kHz;% sorted!
load('./r1rho2and4kHz/R20_WT_15deg_R1r_750_4kHz');R20_750_4kHz = R20_WT_15deg_R1r_750_4kHz;% sorted!

vclist1 = load ('bb_relax_lists_cpmg.txt');
T = 0.05;
vclist = vclist1./T;
numbRes = 59;

vclist_sec_600 = [vclist;4000];%,2000,%./T;
vclist_sec_750 = [vclist;4000];%,2000,%./T;

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

index1(:,1) = 1:1:size(Data600.data,1); % used to unsort later
% sorting:
Peaks20_600_sorted = sortrows(horzcat (resNumbVec_600(:,1),Data600_max.data(:,3),index1));
Peaks20_750_sorted = sortrows(horzcat (resNumbVec_750(:,1),Data750_max.data(:,3),index1));

RelInt_600_sorted = sortrows(horzcat (resNumbVec_600(:,1),Data600.data));
RelInt_750_sorted = sortrows(horzcat (resNumbVec_750(:,1),Data750.data));

% create R1rho peaks 
% R1rhoPeaks_600_2kHz(:,1) = exp(R20_600_2kHz(:,1).*T*(-1)).*Peaks20_600_sorted(:,2)./RelInt_600_sorted(:,21);
% R1rhoPeaks_600_2kHz(:,2) = R1rhoPeaks_600_2kHz(:,1).*R20_600_2kHz(:,2).*T;
% 
% R1rhoPeaks_750_2kHz(:,1) = exp(R20_750_2kHz(:,1).*T*(-1)).*Peaks20_750_sorted(:,2)./RelInt_750_sorted(:,21);
% R1rhoPeaks_750_2kHz(:,2) = R1rhoPeaks_750_2kHz(:,1).*R20_750_2kHz(:,2).*T;

R1rhoPeaks_600_4kHz(:,1) = exp(R20_600_4kHz(:,1).*T*(-1)).*Peaks20_600_sorted(:,2)./RelInt_600_sorted(:,21);
R1rhoPeaks_600_4kHz(:,2) = R1rhoPeaks_600_4kHz(:,1).*R20_600_4kHz(:,2).*T;

R1rhoPeaks_750_4kHz(:,1) = exp(R20_750_4kHz(:,1).*T*(-1)).*Peaks20_750_sorted(:,2)./RelInt_750_sorted(:,21);
R1rhoPeaks_750_4kHz(:,2) = R1rhoPeaks_750_4kHz(:,1).*R20_750_4kHz(:,2).*T;


%%%%%%%%
% From one list of peaks with err to three lists of peaks:
zeros59=zeros(size(R1rhoPeaks_600_4kHz,1),1);

% R1rhoPeaks_600_2kHz_1=zeros59;
% R1rhoPeaks_600_2kHz_2=zeros59;
% R1rhoPeaks_600_2kHz_3=zeros59;
% R1rhoPeaks_600_2kHz_unsort_1=zeros59;
% R1rhoPeaks_600_2kHz_unsort_2=zeros59;
% R1rhoPeaks_600_2kHz_unsort_3=zeros59;
% R1rhoPeaks_750_2kHz_1=zeros59;
% R1rhoPeaks_750_2kHz_2=zeros59;
% R1rhoPeaks_750_2kHz_3=zeros59;
% R1rhoPeaks_750_2kHz_unsort_1=zeros59;
% R1rhoPeaks_750_2kHz_unsort_2=zeros59;
% R1rhoPeaks_750_2kHz_unsort_3=zeros59;

R1rhoPeaks_600_4kHz_1=zeros59;
R1rhoPeaks_600_4kHz_2=zeros59;
R1rhoPeaks_600_4kHz_3=zeros59;
R1rhoPeaks_600_4kHz_unsort_1=zeros59;
R1rhoPeaks_600_4kHz_unsort_2=zeros59;
R1rhoPeaks_600_4kHz_unsort_3=zeros59;
R1rhoPeaks_750_4kHz_1=zeros59;
R1rhoPeaks_750_4kHz_2=zeros59;
R1rhoPeaks_750_4kHz_3=zeros59;
R1rhoPeaks_750_4kHz_unsort_1=zeros59;
R1rhoPeaks_750_4kHz_unsort_2=zeros59;
R1rhoPeaks_750_4kHz_unsort_3=zeros59;


for i=1:size(R1rhoPeaks_600_4kHz,1)

%     [R1rhoPeaks_600_2kHz_1(i),R1rhoPeaks_600_2kHz_2(i),R1rhoPeaks_600_2kHz_3(i),~] = err2distrib(R1rhoPeaks_600_2kHz(i,1),R1rhoPeaks_600_2kHz(i,2));
%     [R1rhoPeaks_750_2kHz_1(i),R1rhoPeaks_750_2kHz_2(i),R1rhoPeaks_750_2kHz_3(i),~] = err2distrib(R1rhoPeaks_750_2kHz(i,1),R1rhoPeaks_750_2kHz(i,2));
    [R1rhoPeaks_600_4kHz_1(i),R1rhoPeaks_600_4kHz_2(i),R1rhoPeaks_600_4kHz_3(i),~] = err2distrib(R1rhoPeaks_600_4kHz(i,1),R1rhoPeaks_600_4kHz(i,2));
    [R1rhoPeaks_750_4kHz_1(i),R1rhoPeaks_750_4kHz_2(i),R1rhoPeaks_750_4kHz_3(i),~] = err2distrib(R1rhoPeaks_750_4kHz(i,1),R1rhoPeaks_750_4kHz(i,2));

end

for i=1:size(R1rhoPeaks_600_4kHz,1)

    % Peaks20_600_sorted(i,3) contains old unsorted order
%     R1rhoPeaks_600_2kHz_unsort_1(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_2kHz_1(i);
%     R1rhoPeaks_600_2kHz_unsort_2(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_2kHz_2(i);
%     R1rhoPeaks_600_2kHz_unsort_3(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_2kHz_3(i);
% 
%     R1rhoPeaks_750_2kHz_unsort_1(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_2kHz_1(i);
%     R1rhoPeaks_750_2kHz_unsort_2(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_2kHz_2(i);
%     R1rhoPeaks_750_2kHz_unsort_3(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_2kHz_3(i);

    R1rhoPeaks_600_4kHz_unsort_1(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_4kHz_1(i);
    R1rhoPeaks_600_4kHz_unsort_2(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_4kHz_2(i);
    R1rhoPeaks_600_4kHz_unsort_3(Peaks20_600_sorted(i,3)) = R1rhoPeaks_600_4kHz_3(i);

    R1rhoPeaks_750_4kHz_unsort_1(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_4kHz_1(i);
    R1rhoPeaks_750_4kHz_unsort_2(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_4kHz_2(i);
    R1rhoPeaks_750_4kHz_unsort_3(Peaks20_750_sorted(i,3)) = R1rhoPeaks_750_4kHz_3(i);

end



% save R1rhoPeaks_600_2kHz_unsort_1 etc. as ser files using 
% Data600_max, Data750_max as templates

% Data600_R1rho_600_2kHz_1 = Data600_max; Data600_R1rho_600_2kHz_1.data(:,3) = R1rhoPeaks_600_2kHz_unsort_1;
% Data600_R1rho_600_2kHz_2 = Data600_max; Data600_R1rho_600_2kHz_2.data(:,3) = R1rhoPeaks_600_2kHz_unsort_2;
% Data600_R1rho_600_2kHz_3 = Data600_max; Data600_R1rho_600_2kHz_3.data(:,3) = R1rhoPeaks_600_2kHz_unsort_3;

Data600_R1rho_600_4kHz_1 = Data600_max; Data600_R1rho_600_4kHz_1.data(:,3) = R1rhoPeaks_600_4kHz_unsort_1;
Data600_R1rho_600_4kHz_2 = Data600_max; Data600_R1rho_600_4kHz_2.data(:,3) = R1rhoPeaks_600_4kHz_unsort_2;
Data600_R1rho_600_4kHz_3 = Data600_max; Data600_R1rho_600_4kHz_3.data(:,3) = R1rhoPeaks_600_4kHz_unsort_3;
 
% Data750_R1rho_750_2kHz_1 = Data750_max; Data750_R1rho_750_2kHz_1.data(:,3) = R1rhoPeaks_750_2kHz_unsort_1;
% Data750_R1rho_750_2kHz_2 = Data750_max; Data750_R1rho_750_2kHz_2.data(:,3) = R1rhoPeaks_750_2kHz_unsort_2;
% Data750_R1rho_750_2kHz_3 = Data750_max; Data750_R1rho_750_2kHz_3.data(:,3) = R1rhoPeaks_750_2kHz_unsort_3;
 
Data750_R1rho_750_4kHz_1 = Data750_max; Data750_R1rho_750_4kHz_1.data(:,3) = R1rhoPeaks_750_4kHz_unsort_1;
Data750_R1rho_750_4kHz_2 = Data750_max; Data750_R1rho_750_4kHz_2.data(:,3) = R1rhoPeaks_750_4kHz_unsort_2;
Data750_R1rho_750_4kHz_3 = Data750_max; Data750_R1rho_750_4kHz_3.data(:,3) = R1rhoPeaks_750_4kHz_unsort_3;

%%%

% sernameWT_600MHz_2kHz_1 = 'WT_600MHz_peaks31.ser';
% r1rho2ser(Data600_R1rho_600_2kHz_1,sernameWT_600MHz_2kHz_1);
% sernameWT_600MHz_2kHz_2 = 'WT_600MHz_peaks32.ser';
% r1rho2ser(Data600_R1rho_600_2kHz_2,sernameWT_600MHz_2kHz_2);
% sernameWT_600MHz_2kHz_3 = 'WT_600MHz_peaks33.ser';
% r1rho2ser(Data600_R1rho_600_2kHz_3,sernameWT_600MHz_2kHz_3);
% 
% sernameWT_750MHz_2kHz_2kHz_1 = 'WT_750MHz_peaks31.ser';
% r1rho2ser(Data750_R1rho_750_2kHz_1,sernameWT_750MHz_2kHz_2kHz_1);
% sernameWT_750MHz_2kHz_2kHz_2 = 'WT_750MHz_peaks32.ser';
% r1rho2ser(Data750_R1rho_750_2kHz_2,sernameWT_750MHz_2kHz_2kHz_2);
% sernameWT_750MHz_2kHz_2kHz_3 = 'WT_750MHz_peaks33.ser';
% r1rho2ser(Data750_R1rho_750_2kHz_3,sernameWT_750MHz_2kHz_2kHz_3);

%%%

sernameWT_600MHz_4kHz_1 = 'WT_15deg_600MHz_peaks34.ser';
r1rho2ser(Data600_R1rho_600_4kHz_1,sernameWT_600MHz_4kHz_1);
sernameWT_600MHz_4kHz_2 = 'WT_15deg_600MHz_peaks35.ser';
r1rho2ser(Data600_R1rho_600_4kHz_2,sernameWT_600MHz_4kHz_2);
sernameWT_600MHz_4kHz_3 = 'WT_15deg_600MHz_peaks36.ser';
r1rho2ser(Data600_R1rho_600_4kHz_3,sernameWT_600MHz_4kHz_3);

sernameWT_750MHz_4kHz_4kHz_1 = 'WT_15deg_750MHz_peaks34.ser';
r1rho2ser(Data750_R1rho_750_4kHz_1,sernameWT_750MHz_4kHz_4kHz_1);
sernameWT_750MHz_4kHz_4kHz_2 = 'WT_15deg_750MHz_peaks35.ser';
r1rho2ser(Data750_R1rho_750_4kHz_2,sernameWT_750MHz_4kHz_4kHz_2);
sernameWT_750MHz_4kHz_4kHz_3 = 'WT_15deg_750MHz_peaks36.ser';
r1rho2ser(Data750_R1rho_750_4kHz_3,sernameWT_750MHz_4kHz_4kHz_3);

%%%

function [m] = r1rho2ser(in,sername)

m.Assignment = in.textdata(2:60)';
m.w1 = in.data(:,1);
m.w2 = in.data(:,2);
m.Height = in.data(:,3);

fid=fopen(sername,'w');                   % open a new file for writing
fprintf(fid,"Assignment    w1         w2      Height\n");
for i=1:length(m.Assignment)
    fprintf(fid,'%s %4.3f %6.3f +%e\n',char(m.Assignment{i}),m.w1(i),m.w2(i),m.Height(i) );                  % then the data
end
fid=fclose(fid);                            % close the new file..

end

function [x1,y1,z1,avgsol1] = err2distrib(avg1,err1)

syms x y z 
eqns = [x+y+z == avg1*3, (avg1-x)^2+(avg1-y)^2+(avg1-z)^2 == err1^2*2*sqrt(3)];
vars = [x y z];
[solx, soly, solz] = solve(eqns,vars,'Real',true);

xx = double(solx);
yy = double(soly);
zz = double(solz);

x1=xx(1);
y1=yy(1);
z1=zz(1);

avgsol1 = (x1+y1+z1)/3;

end