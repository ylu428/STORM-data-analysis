function Clustering_ratio(minSMLs_1, sample, SMLs_result, parent_path, sample_name)
    count = [];
    ratio = [];
    x = [];
    for i = 1:length(minSMLs_1)
        for j = 1:length(sample)
            
            temp = [sum(SMLs_result{j+1,i+1}(:,2)==1), sum(SMLs_result{j+1,i+1}(:,2)==2), ...
                sum(SMLs_result{j+1,i+1}(:,2)>2), sum(SMLs_result{j+1,i+1}(:,2)==0 & SMLs_result{j+1,i+1}(:,1)>20)];
            count = cat(1,count, temp);
            ratio = cat(1,ratio, temp/sum(temp));
            x = cat(1, x, strcat((SMLs_result(j+1, 1)), ' cluster > ', string(SMLs_result(1, i+1))));
        end
    end
    bar_chart = bar(ratio, 'stacked');
    xticklabels(x);
    ylabel('Ratio');
    xtickangle(30);
    legend('Cluster = 1','Cluster = 2', 'Cluster > 2', 'Cluster = 0','Location', 'southwest')
    saveas(gcf, strcat(parent_path, sample_name, '_clus_result.png'));
    close(gcf)
end