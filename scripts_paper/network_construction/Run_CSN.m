expression_tumor_fileName = 'JDS_lung_expression_2_75_var_only_names.csv';
[tumor,~,name_tumor]=importdata(expression_tumor_fileName); 
tumor_data=tumor.data; 

for i=1:size(tumor_data,2)
    C = csnet(tumor_data, i);
    cand=full(C{1,i});
    writematrix(cand, sprintf('H:\\1_PHD\\CSN-master\\CSN_lung_2_75_sample_%.f.csv', i))
end 
