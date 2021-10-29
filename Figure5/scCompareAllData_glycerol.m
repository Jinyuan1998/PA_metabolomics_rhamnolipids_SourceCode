close all; clear

addpath('../matlabfunction');

dir1 = 'growthRate_glycerol';
fileformat = fullfile(dir1, '*.xlsx');
growth_files = dir(fileformat);
growth_files = {growth_files.name};
growthRateTabl=[];
%
for i=1:length(growth_files)
    T = readtable([dir1 '/' growth_files{i}]);
    a = growth_files{i};
    a = strsplit(a, '_');
    a = a{5};
    T.date = repmat({a}, size(T,1),1);
    if i>1
    growthRateTabl = [growthRateTabl; T];
    else
        growthRateTabl = T;
    end
    clear T
end

c = growthRateTabl.condition;
c = cellfun(@(X) strsplit(X, '-'), c, 'UniformOutput', false);
c1 = cellfun(@(X) X{2}, c, 'UniformOutput', false);
c2 = cellfun(@(X) X{4}, c, 'UniformOutput', false);
growthRateTabl.promoter=c1;
growthRateTabl.Conc = c2;
growthRateTabl.Conc2 = zeros(height(growthRateTabl),1);
growthRateTabl.Conc2(ismember(growthRateTabl.Conc, '10mM'))=10;
growthRateTabl.Conc2(ismember(growthRateTabl.Conc, '20mM'))=20;
growthRateTabl.Conc2(ismember(growthRateTabl.Conc, '50mM'))=50;
growthRateTabl.Conc2(ismember(growthRateTabl.Conc, '100mM'))=100;
clear c c1 c2

%
growthRateT_H2O2 = growthRateTabl(contains(growthRateTabl.condition, 'PA14-PrhlAB-H2O2') | ...
                        contains(growthRateTabl.condition, 'PA14-con-H2O2'), :);
d = unique(  growthRateT_H2O2.date);
con = sort(unique(  growthRateT_H2O2.Conc2));
growth_ratio_tbl_PA14 = [];
for i=1:length(d)
    for j=1:length(con)
    r1_phaseI =  growthRateT_H2O2.phaseI(ismember( growthRateT_H2O2.date, d{i}) & ismember(growthRateT_H2O2.promoter, 'con') & growthRateT_H2O2.Conc2 == con(j));
    r2_phaseI =  growthRateT_H2O2.phaseI(ismember( growthRateT_H2O2.date, d{i}) & ismember(growthRateT_H2O2.promoter, 'PrhlAB') & growthRateT_H2O2.Conc2 == con(j));
    r21_phaseI = r2_phaseI ./ r1_phaseI;
    r1_phaseII =  growthRateT_H2O2.phaseII(ismember( growthRateT_H2O2.date, d{i}) & ismember(growthRateT_H2O2.promoter, 'con') & growthRateT_H2O2.Conc2 == con(j));
    r2_phaseII =  growthRateT_H2O2.phaseII(ismember( growthRateT_H2O2.date, d{i}) & ismember(growthRateT_H2O2.promoter, 'PrhlAB') & growthRateT_H2O2.Conc2 == con(j));
    r21_phaseII = r2_phaseII ./ r1_phaseII;
    x = repmat(con(j), size(r21_phaseI,1),1);
    r = [r21_phaseI r21_phaseII x];
    r = array2table(r, 'VariableNames', {'phaseI' 'phaseII' 'H2O2conc'});
    r.date = repmat(d(i), size(r21_phaseI,1),1);
    growth_ratio_tbl_PA14 =[growth_ratio_tbl_PA14; r];
    end
end


%% glm

%% H2O2 conc as number for ratio                
mdl_H2O2_growth_ratio_phaseI =fitglme(growth_ratio_tbl_PA14, 'phaseI~1+ H2O2conc+(1 | date)','link', 'log')

mdl_H2O2_growth_ratio_phaseII=fitglme(growth_ratio_tbl_PA14, 'phaseII~1+ H2O2conc+(1 | date)','link', 'log') 

%%                
% promoter
dir2 = 'PromoterActivity_glycerol';
fileformat = fullfile(dir2, '*.xlsx');
prmt_files = dir(fileformat);
prmt_files = {prmt_files.name};
prmtActTabl=[];

for i=1:length(prmt_files)
    T = readtable([dir2 '/' prmt_files{i}]);
    a = prmt_files{i};
    a = strsplit(a, '_');
    a = a{5};
    T.date = repmat({a}, size(T,1),1);
    if i>1
    prmtActTabl = [prmtActTabl; T];
    else
        prmtActTabl = T;
    end
    clear T
end

c = prmtActTabl.condition;
c = cellfun(@(X) strsplit(X, '-'), c, 'UniformOutput', false);
c1 = cellfun(@(X) X{2}, c, 'UniformOutput', false);
c2 = cellfun(@(X) X{4}, c, 'UniformOutput', false);
prmtActTabl.promoter=c1;
prmtActTabl.Conc = c2;
prmtActTabl.Conc2 = zeros(height(prmtActTabl),1);
prmtActTabl.Conc2(ismember(prmtActTabl.Conc, '10mM'))=10;
prmtActTabl.Conc2(ismember(prmtActTabl.Conc, '20mM'))=20;
prmtActTabl.Conc2(ismember(prmtActTabl.Conc, '50mM'))=50;
prmtActTabl.Conc2(ismember(prmtActTabl.Conc, '100mM'))=100;
clear c c1 c2
%
prmtActT_H2O2 = prmtActTabl(contains(prmtActTabl.condition, 'PA14-PrhlAB-H2O2') | ...
                        contains(prmtActTabl.condition, 'PA14-con-H2O2'), :);
writetable(prmtActT_H2O2, 'Glycerol_prmtActT_H2O2.xlsx');     

%%

d = unique( prmtActT_H2O2.date);
con = sort(unique( prmtActT_H2O2.Conc2));
prmt_ratio_tbl_PA14 = [];
for i=1:length(d)
    for j=1:length(con)
    r1_phaseI =  prmtActT_H2O2.phaseI_median(ismember(prmtActT_H2O2.date, d{i}) & ismember(prmtActT_H2O2.promoter, 'con') & prmtActT_H2O2.Conc2 == con(j));
    r2_phaseI =  prmtActT_H2O2.phaseI_median(ismember( prmtActT_H2O2.date, d{i}) & ismember(prmtActT_H2O2.promoter, 'PrhlAB') & prmtActT_H2O2.Conc2 == con(j));
    r21_phaseI = r2_phaseI ./ r1_phaseI;
    r1_phaseII = prmtActT_H2O2.phaseII_median(ismember(prmtActT_H2O2.date, d{i}) & ismember(prmtActT_H2O2.promoter, 'con') & prmtActT_H2O2.Conc2 == con(j));
    r2_phaseII = prmtActT_H2O2.phaseII_median(ismember(prmtActT_H2O2.date, d{i}) & ismember(prmtActT_H2O2.promoter, 'PrhlAB') & prmtActT_H2O2.Conc2 == con(j));
    r21_phaseII = r2_phaseII ./ r1_phaseII;
    x = repmat(con(j), size(r21_phaseI,1),1);
    r = [r21_phaseI r21_phaseII x];
    r = array2table(r, 'VariableNames', {'phaseI' 'phaseII' 'H2O2conc'});
    r.date = repmat(d(i), size(r21_phaseI,1),1);
    prmt_ratio_tbl_PA14 =[prmt_ratio_tbl_PA14; r];
    end
end

prmt_ratio_tbl_PA14.strain = repmat({'PA14'}, height(prmt_ratio_tbl_PA14), 1);

%%
%% H2O2 conc as number for ratio                
mdl_H2O2_pr_ratio_phaseI =fitglme(prmt_ratio_tbl_PA14, 'phaseI~1+ H2O2conc+(1 | date)', 'link','log')

mdl_H2O2_pr_ratio_phaseII=fitglme(prmt_ratio_tbl_PA14, 'phaseII~1+ H2O2conc+(1 | date)', 'link', 'log') 

%%
%%
growthRateT_rhlA_H2O2 = growthRateTabl(contains(growthRateTabl.condition, 'rhlA-PrhlAB-H2O2') | ...
                        contains(growthRateTabl.condition, 'rhlA-con-H2O2'), :);
d = unique(  growthRateT_rhlA_H2O2.date);
con = sort(unique(  growthRateT_rhlA_H2O2.Conc2));
growth_ratio_tbl_rhlA = [];
for i=1:length(d)
    for j=1:length(con)
        r1_phaseI =  growthRateT_rhlA_H2O2.phaseI(ismember( growthRateT_rhlA_H2O2.date, d{i}) & ismember(growthRateT_rhlA_H2O2.promoter, 'con') & growthRateT_rhlA_H2O2.Conc2 == con(j));
        r2_phaseI =  growthRateT_rhlA_H2O2.phaseI(ismember( growthRateT_rhlA_H2O2.date, d{i}) & ismember(growthRateT_rhlA_H2O2.promoter, 'PrhlAB') & growthRateT_rhlA_H2O2.Conc2 == con(j));
        if ~isempty(r1_phaseI) &  ~isempty(r1_phaseI)
            r21_phaseI = r2_phaseI ./ r1_phaseI;
            r1_phaseII =  growthRateT_rhlA_H2O2.phaseII(ismember( growthRateT_rhlA_H2O2.date, d{i}) & ismember(growthRateT_rhlA_H2O2.promoter, 'con') & growthRateT_rhlA_H2O2.Conc2 == con(j));
            r2_phaseII =  growthRateT_rhlA_H2O2.phaseII(ismember( growthRateT_rhlA_H2O2.date, d{i}) & ismember(growthRateT_rhlA_H2O2.promoter, 'PrhlAB') & growthRateT_rhlA_H2O2.Conc2 == con(j));
            r21_phaseII = r2_phaseII ./ r1_phaseII;
            x = repmat(con(j), size(r21_phaseI,1),1);
            r = [r21_phaseI r21_phaseII x];
            r = array2table(r, 'VariableNames', {'phaseI' 'phaseII' 'H2O2conc'});
            r.date = repmat(d(i), size(r21_phaseI,1),1);
            growth_ratio_tbl_rhlA =[growth_ratio_tbl_rhlA; r];
        end
    end
end

%% glm
%% H2O2 conc as number for ratio                
mdl_H2O2_growth_ratio_phaseI_rhlA =fitglme(growth_ratio_tbl_rhlA, 'phaseI~1+ H2O2conc+(1 | date)','link', 'log')

mdl_H2O2_growth_ratio_phaseII_rhlA=fitglme(growth_ratio_tbl_rhlA, 'phaseII~1+ H2O2conc+(1 | date)','link', 'log') 

%%                
%
prmtActT_H2O2_rhlA = prmtActTabl(contains(prmtActTabl.condition, 'rhlA-PrhlAB-H2O2') | ...
                        contains(prmtActTabl.condition, 'rhlA-con-H2O2'), :);
                    
d = unique( prmtActT_H2O2_rhlA.date);
con = sort(unique( prmtActT_H2O2_rhlA.Conc2));
prmt_ratio_tbl_rhlA = [];
for i=1:length(d)
    for j=1:length(con)    
        r1_phaseI =  prmtActT_H2O2_rhlA.phaseI_median(ismember(prmtActT_H2O2_rhlA.date, d{i}) & ismember(prmtActT_H2O2_rhlA.promoter, 'con') & prmtActT_H2O2_rhlA.Conc2 == con(j));
        r2_phaseI =  prmtActT_H2O2_rhlA.phaseI_median(ismember( prmtActT_H2O2_rhlA.date, d{i}) & ismember(prmtActT_H2O2_rhlA.promoter, 'PrhlAB') & prmtActT_H2O2_rhlA.Conc2 == con(j));
        if ~isempty(r1_phaseI) & ~isempty(r2_phaseI)
            r21_phaseI = r2_phaseI ./ r1_phaseI;
            r1_phaseII = prmtActT_H2O2_rhlA.phaseII_median(ismember(prmtActT_H2O2_rhlA.date, d{i}) & ismember(prmtActT_H2O2_rhlA.promoter, 'con') & prmtActT_H2O2_rhlA.Conc2 == con(j));
            r2_phaseII = prmtActT_H2O2_rhlA.phaseII_median(ismember(prmtActT_H2O2_rhlA.date, d{i}) & ismember(prmtActT_H2O2_rhlA.promoter, 'PrhlAB') & prmtActT_H2O2_rhlA.Conc2 == con(j));
            r21_phaseII = r2_phaseII ./ r1_phaseII;
            x = repmat(con(j), size(r21_phaseI,1),1);
            r = [r21_phaseI r21_phaseII x];
            r = array2table(r, 'VariableNames', {'phaseI' 'phaseII' 'H2O2conc'});
            r.date = repmat(d(i), size(r21_phaseI,1),1);
            prmt_ratio_tbl_rhlA =[prmt_ratio_tbl_rhlA; r];
        end
    end
end

prmt_ratio_tbl_rhlA.strain = repmat({'rhlA'}, height(prmt_ratio_tbl_rhlA),1);

%%
%% H2O2 conc as number for ratio                
mdl_H2O2_pr_ratio_phaseI_rhlA =fitglme(prmt_ratio_tbl_rhlA, 'phaseI~1+ H2O2conc+(1 | date)', 'link','log')

mdl_H2O2_pr_ratio_phaseII_rhlA=fitglme(prmt_ratio_tbl_rhlA, 'phaseII~1+ H2O2conc+(1 | date)', 'link', 'log') 


%% compare promoter acitivity between wt and -rhlA when no H2O2 was added to see if their baseline are different
T = [prmt_ratio_tbl_PA14;prmt_ratio_tbl_rhlA];
T0 = T(T.H2O2conc == 0, :);
mdl_noH2O2_wt_rhlA = fitglme(T0, 'phaseII ~1+ strain+(1 | date)', 'link', 'log');
[p, h]=ranksum(T0.phaseII(strcmp(T0.strain, 'PA14')), T0.phaseII(strcmp(T0.strain, 'rhlA')));
%%
mdl_alldata_wt_rhlA = fitglme(T, 'phaseII ~1+ strain+H2O2conc+(1 | date)', 'link', 'log')
%%
cmap = [.8 .4 .2;
    .2 .6 .8;
    .7 .2 .4;
    .2 .8 .4];
%%
alpha_poly = 0.3;
l_w = 1.5;
%                    
figure
subplot(1,2,1)
mdl = mdl_H2O2_growth_ratio_phaseII;
tbl = growth_ratio_tbl_PA14;
a0=mdl.Coefficients.Estimate(1);
b=mdl.Coefficients.Estimate(2);

a_lower = mdl.Coefficients.Lower(1);
a_upper = mdl.Coefficients.Upper(1);

b_lower = mdl.Coefficients.Lower(2);
b_upper = mdl.Coefficients.Upper(2);

x_glm = sort(unique(prmtActT_H2O2.Conc2));
y_glm = exp(a0+b*x_glm);
% y_glm = a0+b*x_glm;

y_glm_upper = exp(a_upper+(b_upper)*x_glm);
% y_glm_upper(1) = exp(a_upper);

y_glm_lower = exp(a_lower+(b_lower)*x_glm);

x_poly = [x_glm(1) x_glm(end) x_glm(end) x_glm(1) x_glm(1)];
y_poly = [y_glm_lower(1) y_glm_lower(end) y_glm_upper(end) y_glm_upper(1) y_glm_lower(1)];
lineC = cmap(1,:);
plot(tbl.H2O2conc, tbl.phaseII, '^', 'MarkerEdgeColor', lineC, 'lineWidth', l_w)
hold on
plot(x_glm, y_glm, '-', 'Color', lineC , 'lineWidth', 3)
hold on
p_h = patch(x_poly, y_poly, lineC);
p_h.FaceAlpha = alpha_poly;
p_h.EdgeColor = lineC;
p_h.EdgeAlpha = alpha_poly;
% text(0, 1.8, sprintf('G(PrhlAB)/G(Pcon) = %.2f + %.3f * [H2O2]\n p=%.2f', exp(a0), exp(b), mdl.Coefficients.pValue(2)),'Color', lineC)


mdl = mdl_H2O2_growth_ratio_phaseII_rhlA;
tbl = growth_ratio_tbl_rhlA;
a0=mdl.Coefficients.Estimate(1);
b=mdl.Coefficients.Estimate(2);

a_lower = mdl.Coefficients.Lower(1);
a_upper = mdl.Coefficients.Upper(1);

b_lower = mdl.Coefficients.Lower(2);
b_upper = mdl.Coefficients.Upper(2);

x_glm = sort(unique(prmtActT_H2O2_rhlA.Conc2));
y_glm = exp(a0+b*x_glm);

y_glm_upper = exp(a_upper+(b_upper)*x_glm);

y_glm_lower = exp(a_lower+(b_lower)*x_glm);
x_poly = [x_glm(1) x_glm(end) x_glm(end) x_glm(1) x_glm(1)];
y_poly = [y_glm_lower(1) y_glm_lower(end) y_glm_upper(end) y_glm_upper(1) y_glm_lower(1)];
lineC = cmap(2,:);
plot(tbl.H2O2conc, tbl.phaseII, '^', 'MarkerEdgeColor', lineC, 'lineWidth', l_w)
hold on
plot(x_glm, y_glm, '-', 'Color', lineC , 'lineWidth', 3)
hold on
p_h = patch(x_poly, y_poly, lineC);
p_h.FaceAlpha = alpha_poly;
p_h.EdgeColor = lineC;
p_h.EdgeAlpha = alpha_poly;
set(gca, 'xlim', [-10 110], 'xtick', 0:50:100, ...
            'ylim', [.1 2], 'ytick', 0:1:2)
% set(gca, 'Yscale', 'log') 
% text(0, 0.2, sprintf('G(PrhlAB)/G(Pcon) = %.2f + %.3f * [H2O2]\n p=%.2f', exp(a0), exp(b), mdl.Coefficients.pValue(2)), ...
%         'Color', lineC)
xlabel('[H_2O_2] (mM)', 'fontsize', 12)
ylabel('growth ratio', 'fontsize', 12)
title('Growth in Phase II (Glycerol)', 'fontsize', 14)
legend({'wt', 'wt fit', '','-rhlA', '-rhlA fit', ''}, 'Location', 'Southoutside')

%                    
subplot(1,2,2)

mdl = mdl_H2O2_pr_ratio_phaseII;
tbl = prmt_ratio_tbl_PA14;
a0=mdl.Coefficients.Estimate(1);
b=mdl.Coefficients.Estimate(2);

a_lower = mdl.Coefficients.Lower(1);
a_upper = mdl.Coefficients.Upper(1);

b_lower = mdl.Coefficients.Lower(2);
b_upper = mdl.Coefficients.Upper(2);

x_glm = sort(unique(prmtActT_H2O2.Conc2));
y_glm = exp(a0+b*x_glm);

y_glm_upper = exp(a_upper+(b_upper)*x_glm);

y_glm_lower = exp(a_lower+(b_lower)*x_glm);

x_poly = [x_glm(1) x_glm(end) x_glm(end) x_glm(1) x_glm(1)];
y_poly = [y_glm_lower(1) y_glm_lower(end) y_glm_upper(end) y_glm_upper(1) y_glm_lower(1)];
lineC = cmap(1,:);
plot(tbl.H2O2conc, tbl.phaseII, 'o', 'MarkerEdgeColor', lineC, 'lineWidth', l_w)
hold on
plot(x_glm, y_glm, '-', 'Color', lineC , 'lineWidth', 3)
hold on
p_h = patch(x_poly, y_poly, lineC);
p_h.FaceAlpha = alpha_poly;
p_h.EdgeColor = lineC;
p_h.EdgeAlpha = alpha_poly;
% text(0, 12.5, sprintf('PrhlAB/Pcon = %.2f + %.3f * [H2O2]\n p=%.3f', exp(a0), exp(b), mdl.Coefficients.pValue(2)), 'Color', lineC)
hold on
mdl = mdl_H2O2_pr_ratio_phaseII_rhlA;
tbl = prmt_ratio_tbl_rhlA;
a0=mdl.Coefficients.Estimate(1);
b=mdl.Coefficients.Estimate(2);

a_lower = mdl.Coefficients.Lower(1);
a_upper = mdl.Coefficients.Upper(1);

b_lower = mdl.Coefficients.Lower(2);
b_upper = mdl.Coefficients.Upper(2);

x_glm = sort(unique(prmtActT_H2O2_rhlA.Conc2));
y_glm = exp(a0+b*x_glm);

y_glm_upper = exp(a_upper+(b_upper)*x_glm);

y_glm_lower = exp(a_lower+(b_lower)*x_glm);
x_poly = [x_glm(1) x_glm(end) x_glm(end) x_glm(1) x_glm(1)];
y_poly = [y_glm_lower(1) y_glm_lower(end) y_glm_upper(end) y_glm_upper(1) y_glm_lower(1)];
lineC = cmap(2,:);
plot(tbl.H2O2conc, tbl.phaseII, 'o', 'MarkerEdgeColor', lineC, 'linewidth', l_w)
hold on
plot(x_glm, y_glm, '-', 'Color', lineC , 'lineWidth', 3)
hold on
p_h = patch(x_poly, y_poly, lineC);
p_h.FaceAlpha = alpha_poly;
p_h.EdgeColor = lineC;
p_h.EdgeAlpha = alpha_poly;
set(gca, 'xlim', [-10 110], 'xtick', 0:50:100, ...
            'ylim', [4 13], 'ytick', 0:2:12)
set(gca, 'Yscale', 'log')
% text(0, 4.5, sprintf('PrhlAB/Pcon = %.2f + %.3f * [H2O2]\n p=%.3f', exp(a0), exp(b), mdl.Coefficients.pValue(2)), ...
%         'Color', lineC)
xlabel('[H_2O_2] (mM)', 'fontsize', 12)
ylabel('promoter activity ratio', 'fontsize', 12)
title('Promoter Activity in Phase II (Glycerol)', 'fontsize', 14)  
hold off
legend({'wt', 'wt fit', '', '-rhlA', '-rhlA fit', ''}, 'Location', 'Southoutside')
set(gcf, 'Renderer', 'painters')



%%
prmt_ratio_tbl_PA14.strain = repmat({'PA14'}, height(prmt_ratio_tbl_PA14), 1);
prmt_ratio_tbl_rhlA.strain = repmat({'rhlA'}, height(prmt_ratio_tbl_rhlA), 1);

prmt_ratio_tbl = [prmt_ratio_tbl_PA14; prmt_ratio_tbl_rhlA];

mdl_H2O2_pr_ratio_phaseII_all =fitglme(prmt_ratio_tbl, 'phaseII~ 1 +strain + H2O2conc+(1 | date)', 'link', 'log')
exp(mdl_H2O2_pr_ratio_phaseII_all.Coefficients.Estimate(1) + mdl_H2O2_pr_ratio_phaseII_all.Coefficients.Estimate(3)) - ...
    exp(mdl_H2O2_pr_ratio_phaseII_all.Coefficients.Estimate(1))
%%
growth_ratio_tbl_PA14.strain = repmat({'PA14'}, height(growth_ratio_tbl_PA14), 1);
growth_ratio_tbl_rhlA.strain = repmat({'rhlA'}, height(growth_ratio_tbl_rhlA), 1);
growth_ratio_tbl = [growth_ratio_tbl_PA14; growth_ratio_tbl_rhlA];

mdl_H2O2_grwth_ratio_phaseII_all =fitglme(growth_ratio_tbl, 'phaseII~ 1 +strain + H2O2conc+(1 | date)', 'link', 'log')