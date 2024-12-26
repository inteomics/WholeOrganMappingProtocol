patient_ID_list = main_get_patient_list();
main_fit_SPATIAL_UNNORMALIZED_SFS_multiple_cases(patient_ID_list);
figure(1);clf;main_plot_SPATIAL_UNNORMALIZED_SFS_multiple_cases(patient_ID_list);

%=======================================================================
function [patient_ID_list] = main_get_patient_list()
motherfile              = 'WES_MULTIREGION_BLADDER/';
filename                = 'Gene_list.txt';
fileID                  = fopen([motherfile filename],'r');
T                       = textscan(fileID,'%s','delimiter','\n');
patient_ID_list         = T{1};
end
%====================Fit the spatial, unnormalized SFS of multiple cases
function main_fit_SPATIAL_UNNORMALIZED_SFS_multiple_cases(patient_ID_list)
motherfile      = 'WES_MULTIREGION_BLADDER/Data and results for spatial SFS - unnormalized/';
%-------------------Divides the database into mini-database for each map
mini_database                       = cell(1,2);
mini_database{1}                    = {};
mini_database{2}                    = {};
for i_gene=1:size(patient_ID_list,1)
    real_map_ID                     = extractBetween(patient_ID_list{i_gene,1},1,5);
    if strcmp(real_map_ID{1},'MAP19')
        mini_database{1}{end+1,1}   = patient_ID_list{i_gene,1};
        mini_database{1}{end,2}     = patient_ID_list{i_gene,2};
    elseif strcmp(real_map_ID{1},'MAP24')
        mini_database{2}{end+1,1}   = patient_ID_list{i_gene,1};
        mini_database{2}{end,2}     = patient_ID_list{i_gene,2};
    end
end
%-------------------------------------------------------Fit for each map
table_output_all                        = {};
for i_map=1:2
    map_database                        = mini_database{i_map};
    if i_map==1
        file_ID_fitting_c           = fopen([motherfile 'MAP_19_fitting_process_for_c.txt'],'w');
    elseif i_map==2
        file_ID_fitting_c           = fopen([motherfile 'MAP_24_fitting_process_for_c.txt'],'w');
    end
    %       Find the best parameter c for this map
    [para_c,err]                        = main_fit_SPATIAL_UNNORMALIZED_SFS_c(map_database,file_ID_fitting_c);
    %       Find the parameters a and b for each gene
    table_para_a_and_b                  = cell(1,size(map_database,1));
    for i=1:size(map_database,1)
        fprintf('---Fitting %d/%d\n',i,size(map_database,1));
        patient_ID              = map_database{i,1};
        if ~exist(strcat(motherfile,patient_ID,'.txt'),'file')
            disp('DID NOT FIND FILE IN MOTHER FOLDER!!!');
            continue;
        end
        [vec_para_a_and_b,err]      = main_fit_SPATIAL_UNNORMALIZED_SFS_a_and_b(database,patient_ID,para_c);
        fprintf('   Error=%f\n\n',err);
        table_output_all{end+1,1}   = patient_ID;
        table_output_all{end,2}     = map_database{i,2};
        table_output_all{end,3}     = vec_para_a_and_b(1);
        table_output_all{end,4}     = vec_para_a_and_b(2);
        table_output_all{end,5}     = para_c;
        filename                    = [motherfile 'parameter_' patient_ID '.txt'];
        fileID                      = fopen(filename,'w');
        fprintf(fileID,'%.12f ',[vec_para_a_and_b(1) vec_para_a_and_b(2) para_c]);
        fclose(fileID);
    end
end
%----------------------------------------------Output all parameter sets
table_output_all                    = [{'Mutation ID','Map area with max SFS','a','b','c'};table_output_all];
filename                                = [motherfile 'parameter_all.xlsx'];
writecell(table_output_all,filename);
end
%=======================FIT ONE UNNORMALIZED SPATIAL SFS OF A GIVEN GENE
%========================================================FOR PARAMETER C
function [para_c,err] = main_fit_SPATIAL_UNNORMALIZED_SFS_c(map_database,file_ID_fitting_c)
%-----------------------------------------------Fit the spatial VAF data
func_fit        = @(x) error_SPATIAL_UNNORMALIZED_SFS_c(x,database,map_database,file_ID_fitting_c);
x0              = [200];
options         = optimset('TolX',1);
[para_c,err]    = fminbnd(func_fit,150,2*10^4,options);
end
%======COMPUTE THE ERROR OF A PARAMETER SET FOR UNNORMALIZED SPATIAL SFS
%=======================================================FOR PARAMETERS C
function output = error_SPATIAL_UNNORMALIZED_SFS_c(para_c,database,map_database,file_ID_fitting_c)
if para_c<0
    output                      = Inf;
else
    motherfile              = 'WES_MULTIREGION_BLADDER/Data and results for spatial SFS - unnormalized/';
    table_para_a_and_b          = cell(1,size(map_database,1));
    vec_err                     = zeros(1,size(map_database,1));
    parfor i=1:size(map_database,1)
        if (mod(i,1000)==0)
            fprintf('---%d/%d\n',i,size(map_database,1));
        end
        patient_ID          = map_database{i,1};
        if ~exist(strcat(motherfile,patient_ID,'.txt'),'file')
            disp('DID NOT FIND FILE IN MOTHER FOLDER!!!');
            continue;
        end
        [vec_para_a_and_b,err]  = main_fit_SPATIAL_UNNORMALIZED_SFS_a_and_b(patient_ID,para_c);
        table_para_a_and_b{i}   = vec_para_a_and_b;
        vec_err(i)              = err;
    end
    output                      = sum(vec_err);
end
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nParameter c=%f----->Error=%f\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n',para_c,output);
fprintf(file_ID_fitting_c,'%f %f\n',para_c,output);
end
%=======================FIT ONE UNNORMALIZED SPATIAL SFS OF A GIVEN GENE
%=================================================FOR PARAMETERS A AND B
function [vec_para,err] = main_fit_SPATIAL_UNNORMALIZED_SFS_a_and_b(patient_ID,para_c)
%-----------------------------------Input the patient's spatial VAF data
motherfile  = 'WES_MULTIREGION_BLADDER/Data and results for spatial SFS - unnormalized/';
filename        = strcat(motherfile,patient_ID,'.txt');
T               = readtable(filename,'Delimiter',' ','ReadVariableNames',false);
T(:,end)        = [];
T               = table2array(T);
spatial_SFS     = T;
%-----------------------------------------------Fit the spatial VAF data
func_fit        = @(x) error_SPATIAL_UNNORMALIZED_SFS_a_and_b([x para_c],spatial_SFS);
x0              = [log(para_c) 1];
options                 = optimset('MaxIter',10000,'MaxFunEvals',10000);
[vec_para,err,exitflag] = fminsearch(func_fit,x0,options);
%----------Minimize the parameter set without appreciable error increase
vec_para_tmp    = vec_para;
err_tmp         = err;
flag_halving    = 0;
while ((err_tmp/err)<1.001) & min(vec_para_tmp)>10^-10
    flag_halving    = 1;
    vec_para_tmp    = vec_para_tmp/2;
    err_tmp         = func_fit(vec_para_tmp);
end
if flag_halving==1
    vec_para        = 2*vec_para_tmp;
    err             = func_fit(vec_para);
end
end
%======COMPUTE THE ERROR OF A PARAMETER SET FOR UNNORMALIZED SPATIAL SFS
%=================================================FOR PARAMETERS A AND B
function output = error_SPATIAL_UNNORMALIZED_SFS_a_and_b(vec_para,real_spatial_SFS)
if min(vec_para)<0
    output          = Inf;
else
    index_end       = length(real_spatial_SFS)-1;
    SFS_EXPECTED    = SFS_SPATIAL_UNNORMALIZED_expected(vec_para,index_end);
    output          = norm(SFS_EXPECTED-real_spatial_SFS,1);
end
end
%==========================COMPUTE THE EXPECTED UNNORMALIZED SPATIAL SFS
function vec_SFS_expected = SFS_SPATIAL_UNNORMALIZED_expected(vec_para,index_end)
%-----------------------------------------------------Get the parameters
a                   = vec_para(1);
b                   = vec_para(2);
c                   = vec_para(3);
%------------------------------------Compute the theoretical spatial SFS
vec_SFS_expected    = zeros(1,index_end+1);
for i=0:index_end
    if i==0
        value       = exp(a) / c;
    else
        value       = exp(a) * (( b/(b+a) )^i) * gamcdf(a+b,i,1) / c;
    end
    vec_SFS_expected(i+1)   = value;
end
end
%==========================Plot the fitted spatial SFS of multiple cases
function main_plot_SPATIAL_UNNORMALIZED_SFS_multiple_cases(patient_ID_list)
motherfile      = 'WES_MULTIREGION_BLADDER/Data and results for spatial SFS - unnormalized/';
for i=length(patient_ID_list):-1:1
    patient_ID  = patient_ID_list{i};
    if (~exist(strcat(motherfile,patient_ID,'.txt'),'file'))||(~exist([motherfile 'parameter_' patient_ID '.txt'],'file'))
        continue;
    end
    fprintf('---Plotting %d/%d\n',i,length(patient_ID_list));
    %-----------------------------------Input the patient's spatial VAF data
    filename        = strcat(motherfile,patient_ID,'.txt');
    T               = readtable(filename,'Delimiter',' ','ReadVariableNames',false);
    T(:,end)        = [];
    T               = table2array(T);
    spatial_SFS     = T;
    plot_SPATIAL_UNNORMALIZED_SFS(patient_ID,spatial_SFS);
end
end
%============================PLOT THE FITTED SPATIAL SFS OF A GIVEN GENE
function plot_SPATIAL_UNNORMALIZED_SFS(patient_ID,spatial_SFS)
clf;
%   Blue
color_blue      = [0 0.4470 0.7410];
%   Orange
color_orange    = [0.8500 0.3250 0.0980];
%   Yellow
color_yellow    = [0.9290 0.6940 0.1250];
%   Purple
color_purple    = [0.4940 0.1840 0.5560];
%   Green
color_green     = [0.4660 0.6740 0.1880];
%   Cyan
color_cyan      = [0.3010 0.7450 0.9330];
%   Red
color_red       = [0.6350 0.0780 0.1840];
%-----------------------------------Input the patient's spatial VAF data
motherfile  = 'WES_MULTIREGION_BLADDER/Data and results for spatial SFS - unnormalized/';
filename        = [motherfile 'parameter_' patient_ID '.txt'];
vec_para        = readtable(filename,'Delimiter',' ','ReadVariableNames',false,'MultipleDelimsAsOne',true);
vec_para        = table2array(vec_para);
%--------------------------------------------Compute the theoretical SFS
index_end       = length(spatial_SFS)-1;
expected_SFS    = SFS_SPATIAL_UNNORMALIZED_expected(vec_para,index_end);
filename        = strcat(motherfile,patient_ID,'_fitted.txt');
filename=strcat('/Users/dinhngockhanh/Downloads/SFS/',patient_ID,'_fitted.txt');
writematrix(expected_SFS,filename,'delimiter',' ');
%-----------------------------------------------------------Plot the SFS
scatter([0:index_end],spatial_SFS,150,'LineWidth',2,'MarkerFaceColor',color_blue,'MarkerEdgeColor',color_blue,'DisplayName','Data');hold on
plot([0:index_end],expected_SFS,'LineWidth',3,'Color',color_orange,'DisplayName','Theory');
%----------------------------------------------------------Plot settings
legend('Location','northeast');
title_text      = [patient_ID ': a=' num2str(vec_para(1)) ', b=' num2str(vec_para(2)) ', c=' num2str(vec_para(3))];
title(title_text);
xlim([0 index_end]);
set(gca,'FontSize',20);
filename        = [motherfile patient_ID '.png'];
saveas(gcf,filename);
end