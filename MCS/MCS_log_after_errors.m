clear all
close all
clc

   
numiter = 1000; 


true_params = [96.6974793601697	2.28792682457840e-07	0.127731095456701...
    0.307434121518132	0.00239697504494893	4.22687281374013...
 	0.00563627165031936	0.0270642511635155	12725.3611292774	0.679667862013611]; %estimated values 
        
 
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:300;


CD4_data_time = round([1.9,18.0,32.1,49.1,90.1,257.1]).*10+1;
CD8_data_time = round([1.9,18.0,32.1,49.0,89.8,255.4]).*10+1;
HIV_data_time  = round([1.7,5.7,9.7,13.6,17.9,20.8,24.9,27.9,32.0,41.1,49.1,67.3,96.1,178.5,259.4]).*10+1;

initial_cond = [918, 0, 583, 1003];



[~, y_trp] = ode23s(@(t,y)HIV_Within_Host(y,true_params),tforward,initial_cond);


Model_CD4_trp = y_trp(CD4_data_time, 1) +  y_trp(CD4_data_time, 2);
Model_CD8_trp = y_trp(CD8_data_time, 3);
Model_HIV_trp = y_trp(HIV_data_time, 4);


 noiselevel = [0, 0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'lambda', 'beta', 'd',  'delta', 'Psi', 'lambda_z', 'b', 'd_z', 'pi', 'c'};


for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  ViralData = log10((noiselev*randn(length(HIV_data_time),1)).*Model_HIV_trp + Model_HIV_trp);
  CD4Data = log10((noiselev*randn(length(CD4_data_time),1)).*Model_CD4_trp + Model_CD4_trp);
  CD8Data = log10((noiselev*randn(length(CD8_data_time),1)).*Model_CD8_trp + Model_CD8_trp);
  
 
            
            lb = [0, 0, 0.01, 0.1, 0, 0, 0, 0.01, 0, 0];
            ub = [100, 1, 1, 1, 10, 100, 1, 1, 5*1e5, 100];
             
            nonlcon = @param_constraints;
      
             k = fmincon(@(k)err_in_data(k,ViralData, CD4Data,CD8Data), true_params,[],[],[],[],lb,ub,nonlcon,optimset('Display','iter',...
             'MaxFunEvals', 4000, 'MaxIter', 4000));
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

end

function error_in_data = err_in_data(k,ViralData, CD4Data,CD8Data)
dt = 0.1;

tforward = 0:dt:300;

CD4_data_time = round([1.9,18.0,32.1,49.1,90.1,257.1]).*10+1;
CD8_data_time = round([1.9,18.0,32.1,49.0,89.8,255.4]).*10+1;
HIV_data_time  = round([1.7,5.7,9.7,13.6,17.9,20.8,24.9,27.9,32.0,41.1,49.1,67.3,96.1,178.5,259.4]).*10+1;


initial_cond = [918, 0, 583, 1003];
 
 [~,y] = ode23s(@(t,y)HIV_Within_Host(y,k),tforward,initial_cond);
 

   if length(y(:,1)) < length(tforward)
        % Handle the case where y_k contains NaN values
         error_in_data = NaN;  % Set the objective to infinity or another suitable value
    else

Model_CD4_f = log10(y(CD4_data_time, 1) +  y(CD4_data_time, 2));
Model_CD8_f = log10(y(CD8_data_time, 3));
Model_HIV_f = log10(y(HIV_data_time, 4));

 


  error_in_data = (1./length(CD4Data))*(sum(((Model_CD4_f - CD4Data)./Model_CD4_f).^2)).^0.5+...
                    (1./length(CD8Data))*(sum(((Model_CD8_f - CD8Data)./Model_CD8_f).^2)).^0.5 +...
                    (1./length(ViralData))*(sum(((Model_HIV_f - ViralData)./Model_HIV_f).^2)).^0.5;  
                        
   end
 end



 function dy = HIV_Within_Host(y,k)

dy = zeros(4,1);

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


dy(1) = lambda - beta*T*V - d*T ;
dy(2) = beta*T*V  - delta*T_i - Psi*T_i*Z;
dy(3) = lambda_z + b*T_i*Z - d_z*Z;
dy(4) = pi*T_i - c*V;

 
 end


function [c,ceq] = param_constraints(k)
c = k(3) - k(4);
ceq = [];

end