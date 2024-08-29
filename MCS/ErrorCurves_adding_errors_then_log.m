
clear all
close all
clc

numiter = 1;

true_params = [96.6974793601697 2.28792682457840e-07 0.127731095456701 0.307434121518132 0.00239697504494893 4.22687281374013...
    0.00563627165031936 0.0270642511635155 12725.3611292774 0.679667862013611];

X = zeros(length(true_params), numiter);

dt = 0.1;
tforward = 0:dt:300;

t_cd4_measure = round([1.9,18.0,32.1,49.1,90.1,257.1]).*10 + 1;
t_cd8_measure = round([1.9,18.0,32.1,49.0,89.8,255.4]).*10 + 1;
t_v_measure = round([1.7,5.7,9.7,13.6,17.9,20.8,24.9,27.9,32.0,41.1,49.1,67.3,96.1,178.5,259.4]).*10 + 1;


initial_cond = [918, 0, 583, 1003];

[~, y_trp] = ode23s(@(t, y) HIV_Within_Host(y, true_params), tforward, initial_cond);

Model_CD4 = y_trp(t_cd4_measure, 1) + y_trp(t_cd4_measure, 2);
Model_CD8 = y_trp(t_cd8_measure, 3);
Model_Viral = y_trp(t_v_measure, 4);

noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE = zeros(length(noiselevel), length(true_params));

for noisei = 1:length(noiselevel)
    rng default
    noiselev = noiselevel(noisei);
    
    ViralData = log10((noiselev * randn(length(t_v_measure), 1)) .* Model_Viral + Model_Viral);
    CD4Data = log10((noiselev * randn(length(t_cd4_measure), 1)) .* Model_CD4 + Model_CD4);
    CD8Data = log10((noiselev * randn(length(t_cd8_measure), 1)) .* Model_CD8 + Model_CD8); % adding errors then log
    
  
    
%     figure(3 * (noisei - 1) + 1)
%     plot(tforward, log10(y_trp(:, 4)), 'k', 'Linewidth', 1.5)
%     hold on
%     plot(tforward(t_v_measure(:)), ViralData, 'ro')
%     hold on
%     ylabel('Viral load','FontSize',14,'FontName','Arial','FontWeight','bold')
%     xlabel('Time (days)','FontSize',14,'FontWeight','bold')
%     title(['\textbf{Noise Level} $\mathbf{\sigma}$: \textbf{' num2str(noiselev) '}'], 'Interpreter', 'latex')
%     set(gca,'FontSize',14,'linewidth',3,'FontWeight','bold')
%     set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',14)

     figure(3 * (noisei - 1) + 1)
    plot(tforward, log10(y_trp(:, 4)), 'k', 'Linewidth', 3)
    hold on
    plot(tforward(t_v_measure(:)), ViralData, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)  
    hold on
    ylabel('Viral load','FontSize',14,'FontName','Arial','FontWeight','bold')
    xlabel('Time (days)','FontSize',14,'FontWeight','bold')
    title(['\textbf{Noise Level} $\mathbf{\sigma} = ' num2str(noiselev*100) '\%$'], 'Interpreter', 'latex')
    set(gca,'FontSize',14,'linewidth',3,'FontWeight','bold')
    set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',14)

    
    figure(3 * (noisei - 1) + 2)
    plot(tforward, log10(y_trp(:, 3)), 'k', 'Linewidth', 3)
    hold on
    plot(tforward(t_cd8_measure(:)), CD8Data, 'ro','MarkerSize', 8, 'LineWidth', 1.5)
    hold on
    ylabel('CD8 cells','FontSize', 14, 'FontName','Arial','FontWeight','bold')
    xlabel('Time (days)','FontSize',14,'FontWeight','bold' )
    title(['\textbf{Noise Level} $\mathbf{\sigma} = ' num2str(noiselev*100) '\%$'], 'Interpreter', 'latex')
    set(gca,'FontSize',14,'linewidth',3,'FontWeight','bold')
    set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',14)

    figure(3 * (noisei - 1) + 3)
    plot(tforward, log10(y_trp(:, 1) + y_trp(:, 2)), 'k', 'Linewidth', 3)
    hold on
    plot(tforward(t_cd4_measure(:)), CD4Data, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5) 
    hold on
    ylabel('CD4 cells','FontSize', 14, 'FontName','Arial','FontWeight','bold')
    xlabel('Time (days)','FontSize',14,'FontWeight','bold' )
    title(['\textbf{Noise Level} $\mathbf{\sigma} = ' num2str(noiselev*100) '\%$'], 'Interpreter', 'latex')
    set(gca,'FontSize',14,'linewidth',3,'FontWeight','bold')
    set(gca, 'YGrid', 'on', 'XGrid', 'on','LineWidth', 2,'fontsize',14)
end

function dy = HIV_Within_Host(y, k)

dy = zeros(4, 1);

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

dy(1) = lambda - beta * T .* V - d * T;
dy(2) = beta * T .* V - delta * T_i - Psi * T_i .* Z;
dy(3) = lambda_z + b * T_i .* Z - d_z * Z;
dy(4) = pi * T_i - c * V;

end
