clc;
clear;
close all;
warning off;

%% Import data

% dirr = 'C:\Users\z5339055\Downloads\Abir_Project\Abir_Project\database';
dirr = '/Volumes/Miad/Abir_Project/database';
addpath(dirr);

File_name = {
    'TCGA-BRCA-rnaexpr.csv'
    'TCGA-LUAD-rnaexpr.csv'
    'TCGA-LUSC-rnaexpr.csv'
    'TCGA-THCA-rnaexpr.csv'
    'TCGA-HNSC-rnaexpr.xlsx'
    'TCGA-LIHC-rnaexpr.xlsx'
    'TCGA-PRAD-rnaexpr.xlsx'
    'TCGA-STAD-rnaexpr.xlsx'
    };

%% Combine all database to concatenate  (8 databases)

data_combine = [];
Y_combine_label = [];

for i = 1:length(File_name)
    
    for PValueCutoff = [0.01 0.04 0.05]
        for FC = [1 2 3 4]
            
            disp(['PValueCutoff: ' num2str(PValueCutoff)]);
            disp(['FC: ' num2str(FC)]);
            
            file_name = File_name{i};
            
            file_directory = [file_name(6:end-12) '_' 'Result_Pvalue_' num2str(PValueCutoff) '_FC_' num2str(FC)];
            mkdir(file_directory)
            
            opts = detectImportOptions(file_name,'NumHeaderLines',0);
            Data = readtable(file_name,opts);
            
            words_row = Data.Properties.VariableNames;
            
            %             words_row = cellfun(@(S) S(1:end-5), words_row, 'Uniform', 0);
            
            id_normal = contains(lower(words_row),'normal');
            id_tumor  = contains(lower(words_row),'tumor');
            
            
            %     words_row(id_normal) = {[file_name(1:end-11) 'normal']};
            %     words_row(id_tumor) = {[file_name(1:end-11) 'tumor']};
            
            words_row(id_normal) = {[file_name(1:end-11) 'normal']};
            words_row(id_tumor) = {[file_name(1:end-11) 'tumor']};
            
            Class_label = words_row(2:end);
            
            data_combine = Data{:,2:end}; %[data_combine Data{:,2:end}];
            Y_combine_label = Class_label'; %[Y_combine_label;Class_label'];
            
            [count,label] = hist(categorical(Y_combine_label),unique(Y_combine_label));
            
            
            %% Data distribution: boxplot
            
            % Quantile normalization
            data_combine(data_combine == 0) = 1e-6;
            
            x_norm = normalization(data_combine','quantile')';
            
            Y_label_boxplot = Y_combine_label;
            
            Y_label_boxplot(contains(Y_label_boxplot,'normal')) = {'normal'};
            label_boxplot = unique(Y_label_boxplot);
            
            x_boxplot = [];
            for ii = 1:length(label_boxplot)
                
                id_label = contains(Y_label_boxplot,label_boxplot{ii});
                x_boxplot = [x_boxplot median(x_norm(:,id_label),2)];
            end
            
            figure(1)
            boxplot(x_boxplot,label_boxplot,...
                'Widths',0.5,'OutlierSize',4,...
                'PlotStyle','traditional','Symbol','ko',...
                'BoxStyle','outline','LabelOrientation','horizontal');
            
            title('Distribution of Median Genes over Quantile Normalized Samples ')
            ylabel('Median');
            set(gca,'FontSize',10,'XTickLabelRotation',45)
            set(gca,'fontname','times')
            
            saveas(gcf,[file_directory '/Fig1.pdf']);
            saveas(gcf,[file_directory '/Fig1.fig']);
            
            %             print('Fig1','-dpdf');
            %% differential expression analysis
            
            % calculation of p-value for combined_data
            
            data_combine(data_combine == 0) = 1e-6;
            
            
            pvalues = [];
            % id_normal = contains(Y_combine_label,'normal');
            
            for j = 1:size(x_norm,1)
                pval = anova1(x_norm(j,:),Y_combine_label,'off');
                pvalues = [pvalues;pval];
            end
            
            pvalue_adj = mafdr(pvalues,'BHFDR',true);
            
            %% t-test
            
            id_normal = contains(Y_combine_label,'normal');
            
            pvalues_ttest = mattest(x_norm(:,id_normal),x_norm(:,~id_normal));
            
            padj = mafdr(pvalues,'BHFDR',true);
            
            pvalues_ttest = padj;
            %%
            % set Fold changing
            FoldChange = FC;
            
            SigStructure = mavolcanoplot(x_norm(:,id_normal),x_norm(:,~id_normal),...
                pvalues_ttest,'Foldchange',FoldChange,'PCutoff',PValueCutoff,'LogTrans',true);
            
            
            features_genes = str2double(SigStructure.GeneLabels);
            
            %%
            % consider the mean
            meanNormal = mean(x_norm(:,id_normal),2);
            meanTumor = mean(x_norm(:,~id_normal),2);
            
            % consider the dispersion
            dispTreated = std(x_norm(:,id_normal),0,2) ./ meanNormal;
            dispUntreated = std(x_norm(:,~id_normal),0,2) ./ meanTumor;
            
            % plot on a log-log scale
            figure(2);
            
            loglog(meanNormal,dispTreated,'r.');
            hold on;
            loglog(meanTumor,dispUntreated,'b.');
            xlabel('log2(Mean)');
            ylabel('log2(Dispersion)');
            legend('Normal','Tumor','Location','southwest');
            ylim([1e-4 1e2])
            
            set(gca,'fontname','times')
            saveas(gcf,[file_directory '/Fig2.pdf']);
            saveas(gcf,[file_directory '/Fig2.fig']);
            
            %%
            % compute the mean and the log2FC
            meanBase = (meanNormal + meanTumor) / 2;
            foldChange = meanNormal ./ meanTumor;
            log2FC = log2(foldChange);
            
            figure(3)
            % plot mean vs. fold change (MA plot)
            mairplot(meanNormal, meanTumor,'Type','MA','Plotonly',true);
            set(get(gca,'Xlabel'),'String','mean of normalized data')
            set(get(gca,'Ylabel'),'String','log2(fold change)')
            
            set(gca,'fontname','times')
            
            saveas(gcf,[file_directory '/Fig3.pdf']);
            saveas(gcf,[file_directory '/Fig3.fig']);
            
            %             print('Fig3','-dpdf');
            %% P value enrichment
            figure(4);
            histogram(pvalues_ttest,100)
            xlabel('P-value')
            ylabel('Frequency')
            title('P-value enrichment')
            
            set(gca,'fontname','times')
            
            saveas(gcf,[file_directory '/Fig4.pdf']);
            saveas(gcf,[file_directory '/Fig4.fig']);
            
            %             print('Fig4','-dpdf');
            %%
            figure(5)
            scatter(log2(mean(x_norm,2)),log2FC,5,pvalues_ttest,'o')
            colormap(flipud(cool(256)))
            colorbar;
            ylabel('log2(Fold Change)')
            xlabel('log2(Mean of normalized data)')
            title('Fold change by FDR')
            ylim([-1.5 2])
            
            set(gca,'fontname','times')
            saveas(gcf,[file_directory '/Fig5.pdf']);
            saveas(gcf,[file_directory '/Fig5.fig']);
            
            %             print('Fig5','-dpdf');
            %% Top features
            
            x_topfeatures_Boxplot = x_norm(features_genes,:);
            
            x_boxplot = [];
            for ii = 1:length(label_boxplot)
                
                id_label = contains(Y_label_boxplot,label_boxplot{ii});
                x_boxplot = [x_boxplot median(x_topfeatures_Boxplot(:,id_label),2)];
            end
            
            figure(6)
            
            boxplot(x_boxplot,label_boxplot,...
                'Widths',0.5,'OutlierSize',4,...
                'PlotStyle','traditional','Symbol','ko',...
                'BoxStyle','outline','LabelOrientation','horizontal');
            
            title('Distribution of Significant Genes over Quantile Normalized Samples ')
            ylabel('Median');
            set(gca,'FontSize',10,'XTickLabelRotation',45)
            set(gca,'fontname','times')
            
            saveas(gcf,[file_directory '/Fig6.pdf']);
            saveas(gcf,[file_directory '/Fig6.fig']);
            
            %             print('Fig6','-dpdf');
            %% tSNE plot
            
            disp('t-SNE ...')
            
            options = statset('MaxIter',1000);
            Y_tsne = tsne(x_norm','NumDimensions',2,'Algorithm',...
                'exact','Distance','euclidean','Options',options); %euclidean;
            
            figure(7)
            clr = ['rb'];
            
            gscatter(Y_tsne(:,1),Y_tsne(:,2),Y_combine_label,clr,'.',10);
            
            xlabel('t-SNE 1')
            ylabel('t-SNE 2')
            
            set(gca,'FontName', 'Times New Roman')
            
            saveas(gcf,[file_directory '/tsne.pdf']);
            saveas(gcf,[file_directory '/tsne.fig']);
            
            disp('done!')
            
            %%
            Result = SigStructure;
            GeneName = Data.Gene_ID(features_genes);
            
            Result.GeneAll = regexprep(GeneName,'.\d+$','');
            
            Result.Gene_UpRegulated = GeneName(SigStructure.FoldChanges > FoldChange);
            Result.Gene_DownRegulated = GeneName(SigStructure.FoldChanges < FoldChange);
            
            struct2csv(Result,[file_directory '/Result_' file_name(6:end-12) '.csv']);
            
            %% BoxPlot of top 10 important features
            idx_top = features_genes(SigStructure.FoldChanges > FoldChange);
            
            if isempty(idx_top)
                disp(['there is no significant gene for FC:' num2str(FC) ' and PValueCutoff:' num2str(PValueCutoff)]);
            else
                
                if length(idx_top)>10
                    idx_top_10 = idx_top(1:10);
                else
                    idx_top_10 = idx_top;
                end
                
                
                xx_scale = normalize(data_combine','range')';
                top_gene = xx_scale(idx_top_10,:);
                
                top_gene_name = Data.Gene_ID(idx_top_10);
                top_gene_name = regexprep(top_gene_name,'.\d+$','');
                
                figure(8)
                
                grouplabel = top_gene_name;
                labels = ["normal" "tumor"];
                
                id_control = contains(Y_combine_label,'normal');
                id_tumor = contains(Y_combine_label,'tumor');
                
                x_group_boxplot = {top_gene(:,id_control)',top_gene(:,id_tumor)'};
                
                boxplotGroup(x_group_boxplot,'Colors','br','PrimaryLabels',labels,'SecondaryLabels',grouplabel,...
                    'interGroupSpace',2,'groupLabelType','Vertical');
                
                ylabel('Gene Expresion');
                
                title(file_name(6:end-12))
                
                set(gca,'FontSize',10,'XTickLabelRotation',45)
                set(gca,'fontname','times')
                saveas(gcf,[file_directory '/topgenesGroupized.pdf']);
                saveas(gcf,[file_directory '/topgenesGroupized.fig']);
                
                figure(9)
                boxplot(top_gene',top_gene_name,...
                    'Widths',0.5,'OutlierSize',4,...
                    'PlotStyle','traditional','Symbol','ko',...
                    'BoxStyle','outline','LabelOrientation','horizontal');
                
                title('Distribution of Top 10 Significant Genes according to pvalue and Fold changing')
                ylabel('Gene Expresion');
                set(gca,'FontSize',10,'XTickLabelRotation',45)
                set(gca,'fontname','times')
                saveas(gcf,[file_directory '/topgenes.pdf']);
                saveas(gcf,[file_directory '/topgenes.fig']);
            end
            
            close all
            
        end
    end
    
end


