%% convert all matrix files to .mat

files_csv = dir('*.csv');

for i = 1 : length(files_csv)
    
    file_in = readmatrix(files_csv(i).name);
       
    [path, name, ext] = fileparts(files_csv(i).name);
    
    file_out = sprintf('%s.mat', name);
    
    save(file_out, 'file_in');
    
end

%% create total mean, NMDA mean, HC mean

files_mat = dir('*.mat');


for i = 1 : length(files_mat)
    
    matrix = load(files_mat(i).name);
    
    if i == 1
        
        all_mat = matrix.file_in;
        
    else
        
        all_mat = cat(3, all_mat, matrix.file_in);
        
    end
    
end

out_mat = round(mean(all_mat,3));


%% create MST and MST "plus" (avgdeg = 10) for all NWs

%files_mat = dir("*.mat");

% import first column of allsubjects.csv 

for i = 1 : length(files_mat)
    
    load(files_mat(i).name);
    
    [mst, mstp10] = backbone_wu(file_in, 10);
    
    % add i+67 for nmda subjects
    mstname = sprintf('%s_mst.mat', allsubjects{i,1});
    mstp10name = sprintf('%s_mstp10.mat', allsubjects{i,1});
    mstcsv = sprintf('%s_mst.csv', allsubjects{i,1});
    mstp10csv = sprintf('%s_mstp10.csv', allsubjects{i,1});
    
    save(mstname, 'mst')
    csvwrite(mstcsv, mst)
    save(mstp10name, 'mstp10')
    csvwrite(mstp10csv, mstp10)
    
end

%% Difference networks

% NOS difference for raw networks (positive indicates more in HC)
diff = hc - nmda;

% NOS difference for MST networks
diff_mst = hc_mst - nmda_mst;
% percentage difference change (i.e. normalize to total NOS)
diff_mst_perc = diff./hc;
% get connections that are not in HC or not in NMDA
[rowNoNMDA, indNoNMDA] = find(diff_mst_perc == 1);
[rowNoHC, indNoHC] = find(diff_mst_perc == -Inf);

%% Plot number of non-overlapping connections for different MST avgdeg values
% hc must be the mean hc network, nmda the mean nmda network


for i = 2 : 50
    
    [hc_mst, hc_mstp] = backbone_wu(hc, i);
    [nmdamst, nmda_mstp] = backbone_wu(nmda, i);
    
    diff_mstp = hc_mstp - nmda_mstp;
    diff_mstp_perc = diff_mstp./hc_mstp;
    
    [rowNoNMDA, indNoNMDA] = find(diff_mstp_perc == 1);
    [rowNoHC, indNoHC] = find(diff_mstp_perc == -Inf);
    
    unique(i) = length(rowNoHC);
    
end
    
%% create thresholded graphs

%files_mat = dir("*.mat");

% import first column of allsubjects.csv 

for i = 1 : length(files_mat)
    
    load(files_mat(i).name);
    
    matnorm = weight_conversion(file_in, 'normalize');
    
    mat01 = threshold_absolute(matnorm, 0.1);
    
    % add i+67 for nmda subjects
    mstname = sprintf('%s_t01.mat', allsubjects1{i+67,1});
    mstcsv = sprintf('%s_t01.csv', allsubjects1{i+67,1});
    
    save(mstname, 'mat01')
    csvwrite(mstcsv, mat01)
    
end

%% Node Strength analysis over thresholds

clear
files_mat = dir("*.mat");

counter2 = 1;
for i = 1 : length(files_mat)
    
    counter = 1;
    
    load(files_mat(i).name);
    
    for tr = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.4 0.55]
    
        mat01 = threshold_proportional(file_in, tr);
        
        str = sum(mat01,1);
        
        if counter == 1
            
            str_all = str;
            
        else
            str_all = [str_all; str];
            
        end
        
        counter = counter + 1;
        
    end
    
    if counter2 == 1
    
        str_3d = str_all;
        
    else
        
        str_3d = cat(3,str_3d, str_all);
        
    end
    
    counter2 = counter2 + 1;
    
end


    
 %% analyse for disconnected nodes

frag = zeros(1,10); 

for x = 1:10
    
    frag(x) = length(find(str_3d(x,:,:)==0));
    
end

%% Area under the curve for strength
% based on str_3d matrix

auc_mat = zeros(61,84);

% iterate through subjects
for i = 1:61
    % iterate through nodes
    for x = 1:84
        
        % new strength value is AUC of all thresholds
        auc_mat(i,x) = round(trapz(str_3d(:,x,i)));
        
    end
    
end

%% global clustering coefficient for AUC
% iterate over thresholds, get within-subject mean for nodes

%clear
files_mat = dir("*.mat");

counter2 = 1;
for i = 1 : length(files_mat)
    
    counter = 1;
    
    load(files_mat(i).name);
    
    for tr = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.4 0.55]
    
        mat = threshold_proportional(file_in, tr);
        
        mat = weight_conversion(mat, 'normalize');
        
        cc = clustering_coef_wu(mat);
        
        
        if counter == 1
            
            cc_all = cc;
            
        else
            cc_all = [cc_all, cc];
            
        end
        
        counter = counter + 1;
        
    end
    
    cc_sub = mean(cc_all,2);
    
    if counter2 == 1
    
        cc_compl = cc_sub;
        
    else
        
        cc_compl = [cc_compl, cc_sub];
        
    end
    
    counter2 = counter2 + 1;
    
end

mean(cc_compl(:))
std(cc_compl(:))


%% Variant to get global CC for all thresholds

clear
files_mat = dir("*.mat");

counter2 = 1;
for i = 1 : length(files_mat)
    
    counter = 1;
    
    load(files_mat(i).name);
    
    for tr = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.4 0.55]
    
        mat = threshold_proportional(file_in, tr);
        
        mat = weight_conversion(mat, 'normalize');
        
        cc = clustering_coef_wu(mat);
        
        
        if counter == 1
            
            cc_all = cc;
            
        else
            cc_all = [cc_all, cc];
            
        end
        
        counter = counter + 1;
        
    end
    
    if counter2 == 1
    
        cc_3d = cc_all;
        
    else
        
        cc_3d = cat(3,cc_3d, cc_all);
        
    end
    
    counter2 = counter2 + 1;
    
end

sprintf('Mean and std for all')
mean(cc_3d(:))
std(cc_3d(:))

sprintf('Mean and std for 0.1')
mean(mean(cc_3d(:,1,:)))
std(mean(cc_3d(:,1,:)))

sprintf('Mean and std for 0.55')
mean(mean(cc_3d(:,10,:)))
std(mean(cc_3d(:,10,:)))

%% Hub connectivity
% HC hubs (top 25%, >= 3/4 hub score)
% 4: L.SPG(28) R.SPG(77)
% 3: L.PoCG(21) L.PrCG(23) R.PrCG(72) L.RMFG(26) R.RMFG(75) L.SFG(27) R.SFG(76) L.SMG(30) 
% R.SMG(79) L.IN(34) R.IPG(56) R.LOG(59) R.PCU(73)

hchub = [ 21,23,26,27,28,30,34,56,59,72,73,75,76,77,79 ];

% NMDA hubs
% 4: L.SPG(28) R.SPG(77) R.PCU(73) 
% 3: L.PoCG(21) L.PrCG(23) R.PrCG(72) L.PCU(24) L.RMFG(26) R.RMFG(75) L.SFG(27)
% R.SFG(76) L.IN(34) R.IN(83) R.FG(55) R.IPG(56) R.LOG(59) R.SMG(79)

nmdahub = [ 21,23,24,26,27,28,34,55,56,59,72,73,75,76,77,79,83 ];

% create binary matrix
hubmat = zeros(84,84);
restmat = zeros(84,84);

% set nodes and their rows/columns to 1
for i = hchub
    
    hubmat(i,:) = 1;
    hubmat(:,i) = 1;
    
end

% multiply to only keep hub connections
%hubs = file_in .* hubmat;

% this gives ALL hub connections (hub + feeder + peripheral)
% next step: two mats, one with only hub-hub one with hub-rest
nodes = 1:84;
nonhub = setdiff(nodes,hchub);

for i = nonhub
    
    restmat(i,:) = 1;
    restmat(:,i) = 1;
    
end

% matrix for feeder connections (hub to peripheral)
%feedper = hubs .* restmat;
% matrix for rich-club connections
hubcon = ~restmat;
%richclub = hubs .* hubcon;

%% Hub connectivity AUC

files_mat = dir("*.mat");

for i = 1 : length(files_mat)
    
    load(files_mat(i).name);
    
    c = 1;
    for tr = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.4 0.55]
    
        mat = threshold_proportional(file_in, tr);
        
        hubonly = mat .* hubmat;
        feeder = hubonly .* restmat;
        richclub = hubonly .* hubcon;
        
        all_str(c) = sum(mat(:));
        fp_str(c) = sum(feeder(:));
        rc_str(c) = sum(richclub(:));
        
        fp_rel(c) = fp_str(c) / all_str(c);
        rc_rel(c) = rc_str(c) / all_str(c);
        
        c = c + 1;
        
    end
    
    if i == 1
        
        fp_rel_all = fp_rel;
        rc_rel_all = rc_rel;
        
    else
        
        fp_rel_all = vertcat(fp_rel_all, fp_rel);
        rc_rel_all = vertcat(rc_rel_all, rc_rel);
        
    end
    
end

rc_mean = mean(rc_rel_all,2);
feed_mean = mean(fp_rel_all,2);
periph = 1 - (rc_mean + feed_mean);
all_hub = [rc_mean, feed_mean, periph];

