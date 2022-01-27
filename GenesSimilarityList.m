clc
clear all
%%
File_name = {
    'Result_THCA.csv'
    'Result_STAD.csv'
    'Result_PRAD.csv'
    'Result_LUSC.csv'
    'Result_LUAD.csv'
    'Result_LIHC.csv'
    'Result_HNSC.csv'
    'Result_BRCA.csv'
    };


%%
for i = 1:length(File_name)
    file_name = File_name{i};
    opts = detectImportOptions(file_name,'NumHeaderLines',0);
    Data = readtable(file_name,opts);
   
    Gene_All{i} = Data.GeneAll;
    Gene_UpRegulated{i} = Data.Gene_UpRegulated;
    Gene_DownRegulated{i} = Data.Gene_DownRegulated;
    
end

%% Common List of Pairs
Gene_Regulated = {Gene_All,Gene_UpRegulated, Gene_DownRegulated};

% 1:All means Up and Down 2:Up Regulated 3:Down Regulated
for id_gene = [1 2 3]
    Gene_name = Gene_Regulated{id_gene};
jj = 1;
for j = 1:length(File_name)-1
    for k = j+1:length(File_name)
        Common_list{jj} = Similar_list(Gene_name{j},Gene_name{k});
        Compare_item{jj} = [File_name{j}(8:end-4) '_to_' File_name{k}(8:end-4)];
        jaccard_index = Jaccard(Common_list{jj},Gene_name{j},Gene_name{k});
        
        Comm_Jaccard{jj} = num2str(jaccard_index);
        
        Result(jj).List_Name = Compare_item{jj};
        Result(jj).Jaccard = Comm_Jaccard{jj};
        Result(jj).List = Common_list{jj};
        
       jj = jj + 1; 
    end
end

%%
y_1 = Common_list{1};
for ii = 2:length(Common_list)
    y_0 = Similar_list(y_1,Common_list{ii});
    y_1 = y_0;
end

Common_list{end+1} = y_1;


Result(jj+1).List_Name = ['All'];
Result(jj+1).Jaccard = {''};
Result(jj+1).List = y_1;

T = struct2table(Result);

if id_gene == 1
    writetable(T,'Common_List_All.csv');
elseif id_gene == 2
    writetable(T,'Common_List_Up.csv');
else
     writetable(T,'Common_List_Down.csv');
end

end

%% Required Functions
%% Jaccard Index

function y = Jaccard(common_list,x1,x2)

Intersection  = length(common_list);
Union         =  length(x1)+ length(x2) - Intersection;
Jaccard_index =  Intersection/Union;

y = Jaccard_index; 

end


%% Find similar list
function y = Similar_list(x1,x2)

% a_list = min([length(x1) length(x2)]);

y = [];

if length(x1) == length(x2)
    idx = strcmp(vertcat(x1),vertcat(x2));
    y = [x2(idx)];
else
    
    for j = 1:length(x2)
        if sum(strcmp(x1,x2(j))) > 0
            y =[y;x2(j)];
        end
    end
    
end
end