function p_values=get_p_values(cell1,cell2,test,tail)

if nargin < 4
    tail = 'two';
end

n_cells = numel(cell1);

for i=1:n_cells
    
    switch test
        case 'ttest'
            [~,p]=ttest2(cell1{i},cell2{i});
        case 'ranksum'
            p=ranksum(cell1{i},cell2{i});
        case 'permtest'
            p = myPermTest2_refined(cell1{i},cell2{i},[],'all',1000,tail);
        case 'signrank'
            p=signrank(cell1{i},cell2{i});
    end  
    
    %     p_values{i}=sprintf('p=%.3f',p);
    
    if p<0.001
        p_values{i}='***';
    elseif p<0.01
        p_values{i}='**';
    elseif p<0.05
        p_values{i}='*';
    else
        p_values{i}='n.s.';
    end
    
end

