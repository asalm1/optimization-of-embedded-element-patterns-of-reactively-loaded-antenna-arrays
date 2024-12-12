%%
% MDMB, Plot S-parameters
% ------------------------------------------------------------------------
% 14.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%% Clear
clear
clc
close all

%% Whoere to save results
tikzfile = '.\results\tikz_figures\MDMB_result_sparameters.tex';

%% Import S-parameters
Sobj = sparameters('E:\Koodit\public\data_for_optimization_of_reactive_scatterers\simulation_data\scatterer_ports_terminated\MDMB\MWS-run-0003.s5p');

%% Parse and plot
f = Sobj.Frequencies;

figure

lgnds = {};

for it = 1:5
    sii = squeeze(Sobj.Parameters(it,it,:));
    plot(f*1e-9, 20*log10(abs(sii)), 'LineWidth', 2)
    hold on

    lgnds{end+1} = sprintf('$s_{%d%d}$', it, it);
end

hold off

ylabel('$|s_{ii}|$ (dB)', 'Interpreter','latex')
xlabel('Frequency (GHz)', 'Interpreter', 'latex')

legend(lgnds, 'Interpreter', 'latex')
set(gca, 'FontSize', 11)

grid on

matlab2tikz('filename',tikzfile, 'externalData', true, 'standalone',true);
