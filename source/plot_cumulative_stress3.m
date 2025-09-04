%This code is triggered by the "Calc cumulative stress" menu in the stress_plotter2.
delete('Output_files/*')

%check import:
if isempty(app.coul_dir)
    errordlg('No Coulomb files imported!')
    return
elseif isempty(app.fault_dir)
    errordlg('No fault geometry imported!')
    return
end

%read the coulomb files
coulomb_files = dir(app.coul_dir);
coul_data = {coulomb_files.name};
coul_data = coul_data(3:end);

%prepare the table for cumulative stress
cum_stress = readtable(fullfile(app.coul_dir,coul_data{1}),'Delimiter',',','HeaderLines',2,'ReadVariableNames',false);
cum_stress = cum_stress(:,[1 21]);
cum_stress.Properties.VariableNames = {'id','fault_name'};
cum_stress.coulomb = zeros(length(cum_stress.id),1);

% load interseismic stress
if app.interseismicCB.Value == true
    backslip_data = readtable(fullfile(app.interseis_file),'Delimiter',',','HeaderLines',2,'ReadVariableNames',false);
    annual_stress = (table2array(backslip_data(:,18)))/10; %converted to MPa
    event_dates = [app.startdateSp.Value;app.EventTable.Data.Year;app.enddateSp.Value];
    if app.startdateSp.Value > event_dates(1) || app.enddateSp.Value < event_dates(end)
        errordlg('Incorrect start or end date');
        return
    end
end

%load fault geometry from files
files = ls([fullfile(app.fault_dir) '\*.csv']);
n_faults = size(files,1);
geometry = nan(length(cum_stress.coulomb),15); %[idx,x1-4,y1-4,z1-4,row,col]
idx = nan(n_faults,1);
fault_names = strings(n_faults,1);
plot = true(n_faults,1);

for i = 1:n_faults
    file_name = strrep(files(i,:),' ',''); %extract fault name from csv file names
    fault_name = convertCharsToStrings(file_name(1:end-4));
    f_idx = find(strcmp(fault_name,cum_stress.fault_name));
    if any(strcmp(fault_name,cum_stress.fault_name)) %check if the fault is also in the coulomb stress table
        xyz = readmatrix(fullfile(app.fault_dir,file_name));
        geometry(f_idx(1):f_idx(1)+size(xyz,1)-1,2:15) = xyz;
        geometry(f_idx(1):f_idx(1)+size(xyz,1)-1,1) = i;
        idx(i) = i;
        %app.FaultTable.Data.fault(i) = fault_name;
    end
end
geometry(isnan(geometry(:,2)),:) = [];
app.FaultTable.Data.idx = idx;
app.FaultTable.Data(isnan(app.FaultTable.Data.idx),:) = [];

%get indices of faults in cum_stress matrix
fault_indices = zeros(n_faults,2);
for i = 1:n_faults
    fault_indices(i,1) = find(strcmp(app.FaultTable.Data.fault_names(i),cum_stress.fault_name),1);
    fault_indices(i,2) = find(strcmp(app.FaultTable.Data.fault_names(i),cum_stress.fault_name),1,'last');
end
clearvars idx fault_names plot

% prepare event_table
event_list = app.EventTable.Data;
event_list(~event_list.Plot,:) = []; %remove unticked events from list

% create combination matrix
vectors = {[0 0]}; idx = 1;
for yr = event_list.Year(1):event_list.Year(end)
    sel = find(event_list.Year == yr);
    if ~isempty(sel)
        vectors{idx} = event_list.Scenario(sel)';
        idx = idx + 1;
    end
end
% combine all scenario numbers to matrix: (credit to Luis Mendo @ stackoverflow)
n = numel(vectors);
combs = cell(1,n);
[combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); % the reverse order in these two
combs = cat(n+1, combs{:}); % concat the n n-dim arrays along dimension n+1
combs = reshape(combs,[],n); % reshape to obtain desired matrix

titles = {'year','entire_network'}; %create column heads for the output table
for k = 3:n_faults+2
    titles{k} = app.FaultTable.Data.fault_names{k-2};
end

%% calculate stress
metrics = nan(size(combs,1),5); %stressed_elem, percentage, avg_stress, stress_change, stressdrop
comb_safe = cell(size(combs,1),1);
stress_stats = nan(numel(cum_stress.coulomb),size(combs,1)); %matrix to store stresses of ALL combinations
%initiate plot
exp_fig = false;
stress_plot(app.FaultTable.Data,geometry,cum_stress,app.plot3dCB,yr,i,exp_fig,app.plotdateCB,app.colormapDD)

for i = 1:size(combs,1)
    fprintf('Calculating combination %.0f of %.0f\n',i,size(combs,1));
    combination = combs(i,:); combcounter = 1;
    rupture_stress = nan(length(cum_stress.coulomb),1); %array saving stress of all ruptured elements
    r_idx = 1;
    total_stress = app.prestressSp.Value; %mean stress on entire network
    stress_change = zeros(length(combination),1);
%    stressed_elem = 0; slipped_elem = 1; avg_stress = 1; stress_change = 1; s_prev = 0; s_drop = 0; %pre-allocate metrics vars
    cum_stress.coulomb(:,:) = app.prestressSp.Value; %stress at start
    stress_evo = zeros(app.enddateSp.Value-app.startdateSp.Value+1,n_faults+2); %track stress evolution on faults  
    
    single_event_stats = nan(nnz(app.EventTable.Data.Plot),7); %saves year, mean, min, max stress on rupture plane, and percentage of elements >0, 0.1, 0.5 MPa for each event
    %loop going from start to end date, accumulating stress
    count_events = 0;
    source_CST_total = nan(1000,1);
    for yr = round(app.startdateSp.Value):round(app.enddateSp.Value)
        %coseismic
        current_events = find(event_list.Year == yr);
        if ~isempty(current_events)
            count_events = count_events + 1;
            list_events = find(app.EventTable.Data.Year == yr);
            sel = app.EventTable.Data.Event(list_events);
            event_row = find(strcmp(app.EventTable.Data.Event,sel(combination(combcounter))));
            %fprintf('%s\n',event_list.Event(event_row))
            T = readtable(fullfile(app.coul_dir,coul_data{event_row}),'Delimiter',',','HeaderLines',2,'ReadVariableNames',false);
            T = T(:,[1,5,6,7,8,9,10,11,12,13,14,17,18,19,20,21]);
            T.Properties.VariableNames = {'id','length','strike','dip','latslip','dipslip','sig_right','sig_reverse','normal','coul_right','coul_reverse','spec_rake','spec_coul','element_rake','el_rake_coul','fault_name'};
            %T(1:5,:) %check correct import, debugging only
            T.coul_right = T.coul_right/10; T.coul_reverse = T.coul_reverse/10; T.spec_coul = T.spec_coul/10; T.el_rake_coul = T.el_rake_coul/10; %CONVERT bar TO MPa
            
            dipslip_idx = T.dipslip ~= 0; 
            latslip_idx = T.latslip ~= 0;
            source_idx = find((latslip_idx+dipslip_idx) > 0); %indices of patches that slipped

            %save single event metrics:
            elem_pos = nnz(cum_stress.coulomb(source_idx) > 0)/nnz(cum_stress.coulomb(source_idx))*100; %percentage of positively stressed elements for individual events
            elem_01 = nnz(cum_stress.coulomb(source_idx) > 0.1)/nnz(cum_stress.coulomb(source_idx))*100; %percentage of fault elements > 0.1 MPa
            elem_05 = nnz(cum_stress.coulomb(source_idx) > 0.5)/nnz(cum_stress.coulomb(source_idx))*100; %percentage of fault elements > 0.5 MPa
            single_event_stats(count_events,:) = [yr,mean(cum_stress.coulomb(source_idx)),min(cum_stress.coulomb(source_idx)),max(cum_stress.coulomb(source_idx)),elem_pos,elem_01,elem_05];
            writematrix(single_event_stats,strcat('Output_files/single_event_stats_',num2str(i),'.csv'))

            %save cumulative CST on source fault prior to rupture
            source_CST = cum_stress.coulomb(source_idx);
            source_CST = [yr;source_CST]; %#ok<AGROW>
            source_CST(numel(source_CST)+1:size(source_CST_total,1)) = NaN;
            source_CST_total = [source_CST_total,source_CST]; %#ok<AGROW>
            %writematrix(source_CST,strcat('Output_files/source_fault_CST_',sel(combination(combcounter)),'.csv'))

            rupture_stress(r_idx:r_idx+numel(source_idx)-1) = cum_stress.coulomb(source_idx);
            r_idx = r_idx+numel(source_idx);
            
            %calc metrics for plausibility evaluation:
%            stressed_elem = stressed_elem + nnz(cum_stress.coulomb(source_idx) > 0);
%            slipped_elem_old = slipped_elem;
%            slipped_elem = slipped_elem + numel(source_idx);
%            avg_stress = (slipped_elem_old*avg_stress + numel(source_idx)*mean(cum_stress.coulomb(source_idx)))/slipped_elem;
%            stress_change = (combcounter*stress_change+mean(T.el_rake_coul))/(combcounter+1);
%            s_prev = mean(cum_stress.coulomb(source_idx));

            %accumulate stress
            T.el_rake_coul(source_idx) = 0;
            cum_stress.coulomb = T.el_rake_coul + cum_stress.coulomb;
            % if rb_zero_sd.Value == false %&& sum(T.el_rake_coul) > 0 %exclude empty events (zero stress change)
            %     cum_stress = stress_drop(cum_stress,source_idx,geometry,meandrop_sp); %use bull's eye stress drop fcn
            % elseif rb_zero_sd.Value == true %&& sum(T.el_rake_coul) > 0
            cum_stress.coulomb(source_idx) = 0; %set stress to 0 for slipped segments
            % end
%            s_drop = (slipped_elem_old*s_drop + numel(source_idx)*(s_prev - mean(cum_stress.coulomb(source_idx))))/slipped_elem;
            
            stress_change(combcounter) = mean(cum_stress.coulomb) - total_stress;
            total_stress = mean(cum_stress.coulomb);
            combcounter = combcounter+1;
            
        end
        %interseismic
        if app.interseismicCB.Value == true
            cum_stress.coulomb = cum_stress.coulomb + annual_stress;
        end
        
        for j = 2:n_faults+1
            stress_evo(yr-app.startdateSp.Value+1,j+1) = mean(cum_stress.coulomb(fault_indices(j-1,1)+1:fault_indices(j-1,2)));
        end
        stress_evo(yr-app.startdateSp.Value+1,2) = mean(cum_stress.coulomb);
        stress_evo(yr-app.startdateSp.Value+1,1) = yr;

        %decide what to plot and export:
        switch app.ExpFigDD.Value
            case 'annually'
                exp_fig = true;
                stress_plot(app.FaultTable.Data,geometry,cum_stress,app.plot3dCB,yr,i,exp_fig,app.plotdateCB,app.colormapDD)
            case 'after each event'
                exp_fig = false;
                e = find(app.EventTable.Data.Year == yr);
                if ~isempty(e) && any(app.EventTable.Data.Plot(e))
                    exp_fig = true;
                    stress_plot(app.FaultTable.Data,geometry,cum_stress,app.plot3dCB,yr,i,exp_fig,app.plotdateCB,app.colormapDD)
                end
            case 'end of model'
                if yr == app.enddateSp.Value
                    exp_fig = true;
                    stress_plot(app.FaultTable.Data,geometry,cum_stress,app.plot3dCB,yr,i,exp_fig,app.plotdateCB,app.colormapDD)
                end
            case 'none'
                exp_fig = false;
        end
    end
    
    %EVALUATION METRICS:
    rupture_stress(isnan(rupture_stress)) = [];
    mean_rup_sig = mean(rupture_stress);                        %mean stress on ruptured segments
    sig00 = (nnz(rupture_stress>0)/nnz(rupture_stress)*100);    %percentage of positively stressed patches ruptured
    sig01 = (nnz(rupture_stress>0.1)/nnz(rupture_stress)*100);  %percentage of patches ruptured with stress above 0.1
    sig02 = (nnz(rupture_stress>0.2)/nnz(rupture_stress)*100);
    sig05 = (nnz(rupture_stress>0.5)/nnz(rupture_stress)*100);
    sig1 = (nnz(rupture_stress>1)/nnz(rupture_stress)*100);
    mean_stress_change = mean(stress_change);
    
    %stress_plot(fault_table,geometry,cum_stress,plot3d_cb,vik,yr,i,exp_fig,ANATOLIA)
    fprintf('Combination: %s \n',num2str(combination))
    fprintf('Ruptured segments over 0.1 MPa: %.1f percent\n',sig01)
    fprintf('Mean stress before rupture: %.3f MPa\n',mean_rup_sig)
    fprintf('Mean stress change on fault network: %.3f MPa\n',mean_stress_change)
    %fprintf('Mean stress drop at rupture area: %.3f MPa\n',s_drop)
    fprintf('#########################################################\n')
    
    metrics(i,1:7) = (round([sig00,sig01,sig02,sig05,sig1,mean_rup_sig,mean_stress_change],4));
    comb_safe{i} = num2str(combination);
    % no = strrep(sprintf('%4.0f',i),' ','0');
    % evotable = array2table(stress_evo);
    % evotable.Properties.VariableNames = titles;
    % writetable(evotable,strcat('Output_files\stress_evolution_',no,'.csv'))
    stress_stats(:,i) = cum_stress.coulomb; %save final stress state of each combination
end
stress_plot(app.FaultTable.Data,geometry,cum_stress,app.plot3dCB,yr,i,exp_fig,app.plotdateCB,app.colormapDD)
metrics_data = table(comb_safe,metrics(:,1),metrics(:,2),metrics(:,3),metrics(:,4),metrics(:,5),metrics(:,6),metrics(:,7));
metrics_data.Properties.VariableNames = {'combination','stress00','stress01','stress02','stress05','stress1','mean_stress','stress_change'};
writetable(metrics_data,'Output_files\evaluation_metrics.csv')
clearvars current_events event_row i idx j n n_faults no stress_evo

%save CST on source faults prior to rupture:
source_CST_total(:,1) = [];
writematrix(source_CST_total,'Output_files/source_fault_CST.csv')

%save variables needed for evaluation code
save('Output_files/stress_stats.mat',"stress_stats")
save('Output_files/cum_stress.mat',"cum_stress")
