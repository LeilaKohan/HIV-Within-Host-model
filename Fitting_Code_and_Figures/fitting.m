clear all
close all
clc

global  tforward  ViralData CD4Data immunecell dt



ViralData = [3.552136488,5.061135419,5.973394465,6.326253954,6.044945059,5.558991767,5.203178557,4.885014273,4.697034155,4.619047094,...
    4.504160559,4.407919351,4.384951154,4.167844523,3.876532945]';%log


CD4Data = [905.707196,596.7741935,622.8287841,712.1588089,692.0371701,603.3216963]'; %linear (cells/\mu l)

immunecell = [563.2754342,1207.19603,1132.754342,1088.08933,978.6425711,930.8963959]';%linear (cells/\mu l)




tVLdata = round([1.7,5.7,9.7,13.6,17.9,20.8,24.9,27.9,32.0,41.1,49.1,67.3,96.1,178.5,259.4]');

tCD4data = round([1.9,18.0,32.1,49.1,90.1,257.1]');

timmunecell = round([1.9,18.0,32.1,49.0,89.8,255.4]');


dt = 0.1;
tforward = 0:dt:300;

%params = [lambda beta d delta Psi lambda_z b d_z pi c]

lb = [0,0,0.01,0.1,0,0,0,0.01,0,0];
ub = [100, 1, 1, 1, 10, 100, 1, 1, 5*1e5,100];



%params = [lambda beta d delta Psi lambda_z b d_z pi c T Z V]
%  k=[10,1e-5,0.03,0.24,0.5,10,1e-3,0.03,6112,0.5, 800,500,1000];

%  k = [96.7480338096752 2.50415785876445e-07 0.148291444279105 0.230515411791119...
%      0.00223815724683836 4.17319495303451 0.00371172666899770 0.0121794943463912...
%    12725.3614266140 0.734934380692624 918.243819556931 582.859487488093...
%      1002.60752849915];
%  k=[96.6976282893616	2.34606096885518e-07	0.126274275551617	0.306994563065801	0.00245359801192904	4.22550628300856...
%  	0.00564746576560377	0.0270960915883315	12725.3611243235	0.697382306723785	918.215629296927	582.867318389540...
%     1002.59733811291];
%  k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
%  	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611	918.215708682280	582.867266475109...
%     1002.59735591085];
%  k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
%  	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611	918.215708682280	582.867266475109...
%     1002.59735591085];
%  k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
%  	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611	918.215708682280	582.867266475109...
%     1002.59735591085];

 % k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
 	% 0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611];
 % k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
 	% 0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611];
%  k=[96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
%  	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611]; %estimated values
 k = [96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013...
 	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611];
 %results

 %k = 96.6974793601697	2.28792682457840e-07	0.127731095456701	0.307434121518132	0.00239697504494893	4.22687281374013	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611
initial_cond = [918,0, 583, 1003];



nonlcon = @constraint_function;
 k = fmincon(@(k)err_in_data(k), k,[], [], [], [], lb, ub,nonlcon,optimset('Display','iter'))


[t_r, y_r] = ode23s(@(t,y)HIV_Within_Host(y,k),tforward,initial_cond);


figure(1)
plot(tforward, log10(y_r(:,4)),'b','LineWidth',3)
hold on 
plot(tVLdata, ViralData, 'r.','MarkerSize',25)
set(gca,'FontSize',15,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',15)
title('Viral Load','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Times(days)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Log_{10} V(t)','FontSize',18,'FontName','Arial','FontWeight','bold')



figure(3)
plot(tforward, log10(y_r(:,1)+y_r(:,2)),'b','LineWidth',3)
hold on 
plot(tCD4data, log10(CD4Data), 'r.','MarkerSize',25)
set(gca,'FontSize',15,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',15)
title('Target Cells','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Times(days)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Log_{10} CD4 cells','FontSize',18,'FontName','Arial','FontWeight','bold')



figure(5)
plot(tforward, log10(y_r(:,3)),'b','LineWidth',3)
hold on 
plot(timmunecell, log10(immunecell), 'r.','MarkerSize',25)
set(gca,'FontSize',15,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',15)
title('Immune Cells','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Times(days)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Log_{10} CD8 cells','FontSize',18,'FontName','Arial','FontWeight','bold')

fprintf('lambda = %g\n',  k(1));   
fprintf('beta = %g\n', k(2));
fprintf('d = %g\n', k(3));
fprintf('delta = %g\n',  k(4));
fprintf('Spi = %g\n', k(5));
fprintf('lambda_z = %g\n',  k(6));
fprintf('b = %g\n',  k(7));
fprintf('d_z = %g\n',  k(8));
fprintf('pi = %g\n',  k(9));
fprintf('c = %g\n', k(10));



% Define constraints
function [c, ceq] = constraint_function(k)
    c= k(3)- k(4);
    ceq= [];
    
end



 
function error_in_data = err_in_data(k)

global  tforward  ViralData CD4Data immunecell

initial_cond = [918.215708682280,0, 582.867266475109, 1002.59735591085];

[~,y] = ode23s(@(t,y)HIV_Within_Host(y,k),tforward,initial_cond);

t_v_measure = round([1.7,5.7,9.7,13.6,17.9,20.8,24.9,27.9,32.0,41.1,49.1,67.3,96.1,178.5,259.4]).*10 + 1;

t_cd4_measure = round([1.9,18.0,32.1,49.1,90.1,257.1]).*10 + 1;

t_immune_measure = round([1.9,18.0,32.1,49.0,89.8,255.4]).*10 + 1;





Model_Viral = log10(y(t_v_measure(:), 4));
Model_CD4 = log10(y(t_cd4_measure(:), 1)+y(t_cd4_measure(:), 2));
Model_immune = log10(y(t_immune_measure(:), 3));



% error_in_data = (1./length(ViralData))*(sum(((Model_Viral - ViralData)./Model_Viral).^2))+...
%                  (1./length(CD4Data))*(sum(((Model_CD4 - log10(CD4Data))./Model_CD4).^2))+...
%                   3*(1./length(immunecell))*(sum(((Model_immune - log10(immunecell))./Model_immune).^2));

% 
%  

error_in_data = (1./length(ViralData))*(sum(((Model_Viral - ViralData)./Model_Viral).^2)).^0.5+...
                 (1./length(CD4Data))*(sum(((Model_CD4 - log10(CD4Data))./Model_CD4).^2)).^0.5 +...
                  (1./length(immunecell))*(sum(((Model_immune - log10(immunecell))./Model_immune).^2)).^0.5;


 end




function dy = HIV_Within_Host(y,k)

dy = zeros(4,1);

%params = [lambda beta d delta Psi lambda_z b d_z pi c]

lambda = k(1);
beta = k(2);
d = k(3);
delta = k(4);
Psi = k(5);
lambda_z = k(6);
b = k(7);
d_z = k(8);
pi = k(9);
c = k(10);



T = y(1);
T_i = y(2);
Z = y(3);
V = y(4);


dy(1) = lambda - beta*T.*V - d*T ;
dy(2) = beta*T*V  - delta*T_i - Psi*T_i.*Z;
dy(3) = lambda_z + b*T_i.*Z - d_z*Z;
dy(4) = pi*T_i - c*V;

 
end