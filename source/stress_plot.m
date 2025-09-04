function stress_plot(FaultTable,geometry,cum_stress,plot3dCB,yr,i,exp_fig,plotdateCB,colormapDD)
    load("vik.mat",'vik')
    anatolia = load('anatolia.mat','ANATOLIA');
    figure(3)
    clf(figure(3))
    gcf = figure(3);
    % set(gcf,'WindowState','maximized','Color',[1 1 1])
    set(gcf,'Position', [100, 100, 1200, 600]);
    hold on
    for j = 1:length(FaultTable.idx)
        fault_name = FaultTable.fault_names{j};
        stress_idx = find(strcmp(fault_name,cum_stress.fault_name));
        stress = cum_stress.coulomb(stress_idx(1):stress_idx(end));
        for c = 1:numel(stress_idx)
            x = geometry(c+min(stress_idx)-1,2:5);
            y = geometry(c+min(stress_idx)-1,6:9);
            z = geometry(c+min(stress_idx)-1,10:13);
            if FaultTable.plot(j) == true
                patch(x,y,z,stress(c),'LineStyle','--','EdgeColor','none'); %'EdgeColor','none' 'EdgeColor',[.5 .5 .5]
            end
        end
        surf = find(geometry(stress_idx(1):stress_idx(end),10)==0);
        trace = nan(numel(surf),3);
        for c = 1:numel(surf)
            trace(c,1:3) = geometry(c+stress_idx(1)-1,[2 6 10]);
        end
        trace(c,1:3) = geometry(c+stress_idx(1)-1,[5 9 13]);
        if FaultTable.plot(j) == true
            plot3(trace(:,1),trace(:,2),trace(:,3),'k','LineWidth',2)
        end
    end
    axis('equal')
    if plot3dCB.Value == true
        view(3)
    end
    
    %colorbar_limits = [-cb_limit_sp.Value cb_limit_sp.Value];
    colorbar_limits = [-1 1];
    switch colormapDD.Value
        case 'vik'
            colormap(vik)
        case 'anatolia'
            colormap(anatolia)
            set(gca,'Color',[.8 .8 .8]);
    end
    
    cb = colorbar('northoutside');
    clim(colorbar_limits)
        
    title(cb,'Coulomb stress on faults (MPa)','FontSize',14);
    xlabel('UTM x')
    ylabel('UTM y')
    zlabel('Depth (m)')
    if plotdateCB.Value == true
        annotation(figure(3),'textbox',[.2 .70 .045 .05],'String',num2str(yr),'FontSize',16,'FontWeight','bold','FitBoxToText','off');
    end
    hold off
    if exp_fig == true
        exp_filename = strcat('Output_files\cumulative_stress_',num2str(i),'_',num2str(yr),'.png');
        exportgraphics(gca,exp_filename,'Resolution',400);
    end
end