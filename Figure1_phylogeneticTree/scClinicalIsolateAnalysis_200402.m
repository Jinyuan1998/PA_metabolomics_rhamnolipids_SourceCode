clear; close all;

removeStrain = {'RMI'}; 
myversion = 2021;
%% Load the core SNP (DNA) matrix and Variance co-Variance Matrix information
% read core SNP matrix (core DNA sequences) used to build full phylogenetic Tree
coreSNPtable = readtable('snpMatrixT_31genomes.csv','ReadVariableNames', true, 'ReadRowNames', true); 
refPA14table = readtable('referencePA14.csv','ReadVariableNames', true); 

coreSNPtable.Properties.VariableNames(find(ismember(coreSNPtable.Properties.VariableNames, 'PA14'))) = {'pa_UCBPP_PA14'};
coreSNPtable.Properties.VariableNames(find(ismember(coreSNPtable.Properties.VariableNames, 'W70332'))) = {'W70322'};
sStrain = coreSNPtable.Properties.VariableNames;

snpMatrix = table2array(coreSNPtable);
ZeroIdx = find(sum(snpMatrix, 2)==0);
snpMatrix(ZeroIdx, :) = [];
coreSNPtable(ZeroIdx, :) = [];

snpMatrix = snpMatrix';
columnNames = coreSNPtable.Properties.VariableNames;
rowNames = coreSNPtable.Properties.RowNames;
% get unique gene names from snp names 
gNameinsnp = repmat({''}, length(rowNames), 1);
for i=1:length(rowNames)
    x = strsplit(rowNames{i}, '_');
    gNameinsnp{i} = x{1};    
end
ugNames = unique(gNameinsnp);
% usedgenes = zeros(size(refPA14table, 1), 1);
Alllength = 0;
for i = 1:length(ugNames)
    idx = find(ismember(refPA14table.Var1(:), ugNames{i}));
    if ~isempty(idx)
        Alllength = Alllength + refPA14table.length(i);
    end
end
% use the ratio to correct for phylogenetic scale in phylogenetic tree.
r = size(snpMatrix, 2) / Alllength;

RHtbl=readtable('rhramnolipidsProduction.xlsx','ReadVariableNames', true);
% 
RHtbl.bodysites(contains(RHtbl.strains, 'PA7') | contains(RHtbl.strains, 'PAO1')) = {'other type strain'; 'other type strain'};
% 
RHtbl.bodysites(contains(RHtbl.strains, 'PA14')) = {'wildtype'};
[~, ~, c]=intersect(sStrain, RHtbl.strains);
RHtbl = RHtbl(c, :);

%% creat color map for body sites
RHtbl.bodysites=categorical(RHtbl.bodysites);
categories(RHtbl.bodysites)
RHtbl.bodysites=reordercats(RHtbl.bodysites,{'wildtype','other type strain','Blood', ...
                                            'Body fluid', 'Pubic bone', 'Sputum', ...
                                             'Stool', 'Tissue', 'Urine'})
nUsites = cellstr(nominal(unique(RHtbl.bodysites)));
% bodysitecolormap = [ 204 51 51;                                                           
%                     105 210 210;
%                     200 200 200;
%                     255 153 0; 
%                     204 204 0;
%                     0 102 204;
% %                     0 102 102;  
%                     102 53 201; 
%                     0 114 15;            
%                     51 51 51] / 255;
bodysitecolormap = [ 51 51 51;
                    200 200 200;
                    204 51 51;                      
                    105 210 210;                   
                    255 153 0; 
                    204 204 0;
                    0 102 204;
%                     0 102 102;  
                    102 53 201; 
                    0 114 15 ] / 255;
figure
A =[ones(size(bodysitecolormap,1), 1), (1:size(bodysitecolormap,1))' * 0.1];
for i=1:size(A,1)
plot(A(i,1), A(i,2), 'o', ...
       'MarkerFaceColor', bodysitecolormap(i,:), ...
       'MarkerEdgeColor', bodysitecolormap(i,:), ...
       'MarkerSize',10)
hold on
end
hold off
legend(nUsites)
% create color file for bar
ColorCodeBodySites = zeros(size(RHtbl,1), 3);

% RHtbl.bodysites = categorical(RHtbl.bodysites);
for i = 1:length(nUsites)
    siteIdx = find(RHtbl.bodysites  == nUsites{i});
%     ColorCodeBodySites = [ColorCodeBodySites; repmat(bodysitecolormap(i,:), length(siteIdx), 1)];
    ColorCodeBodySites(siteIdx,:) = repmat(bodysitecolormap(i,:), length(siteIdx), 1);
    clear siteIdx
end
ColorCodeBodySitesT = table(ColorCodeBodySites, 'VariableNames', {'color'});
RHtbl = [RHtbl, ColorCodeBodySitesT];

%% build a phylogenetic tree using core genome and plot phenotype
Dist = pdist(snpMatrix, 'hamming');
phylotree = seqneighjoin(Dist, 'equivar', columnNames);

ind = getbyname(phylotree,'PAO1');
names = get(phylotree, 'LEAFNAMES');
% nameT = array2table(names);
% writetable(nameT, 'namesOrderTree.csv')
t2 = prune(phylotree, ind);
ind = getbyname(t2,'PA7');
t3 = prune(t2, ind);
view(t3)
phytreewrite('coregeneTree', phylotree)
phytreewrite('coregeneTree29', t3)
% phytreewrite('coregeneTreePatricID', phylotree2, 'BranchNames', false)
%view(phylotree)
% figure
hold on
h = plot(phylotree,'Type', 'square');
set(gca, 'XLim', [-.05, 1])
axis('off')
hold on
namesTlb = table(names, 'VariableNames', {'strains'},'Rownames', names);
rhm = categorical(RHtbl.rhInGlycerol);
rhm(rhm == '+')='3';
rhm(rhm == '+/-')='1.5';
rhm(rhm == '-')='0.01';
rhm(isundefined(rhm))='0.1';
if myversion == 2015
    newRhm = str2num(char(rhm));
else
    newRhm = double(string(rhm));
end

for i = 1:length(ColorCodeBodySites)
    ci = find(ismember(RHtbl.strains, names{i}));
    plot(h.LeafDots.XData(i), h.LeafDots.YData(i), 'wo',...
                     'MarkerFaceColor',  RHtbl.color(ci,:), ...
                     'MarkerEdgeColor',  RHtbl.color(ci,:), ...
                     'MarkerSize', 10);
end

xCCord = h.LeafDots.XData +.11;
xBCord = h.LeafDots.XData+.14;
xSCord = h.LeafDots.XData+.17;
rhmSize = newRhm * 4;
mxSCord = max(xSCord)+0.05;
for i = 1:length(ColorCodeBodySites) 
    idx = find(contains(RHtbl.strains, h.leafNodeLabels(i).String));
    plot(mxSCord+0.15, h.LeafDots.YData(i), 'ko',...
                            'MarkerFaceColor', [0.9 0.7 0.7 ],'MarkerEdgeColor', [0.9 0.7 0.7 ], ...
                            'MarkerSize', rhmSize(idx));
    hold on
    box off

end
hold off
box off
axis off

