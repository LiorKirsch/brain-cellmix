%Translate the affy IDs to entrez and gene symbols
% using the output of "David gene conversion web tool" 
%  http://david.abcc.ncifcrf.gov/conversion.jsp

z = textscan( fopen('affy_to_entrez.csv'), '%q %q','delimiter',',','HeaderLines',1);

affy = z{1};
entrez = z{2};
uniq_affy = unique(affy);
uniq_entrez = unique(entrez);
affy2entrez = false(length(uniq_affy), length(uniq_entrez));

for i=1:length(uniq_affy)
    currrent_affy = uniq_affy{i};
    inds_current_affy = strcmp(affy, currrent_affy);
    current_entrez = entrez(inds_current_affy);
    
    [~ ,inds_current_affy] = ismember(current_entrez, uniq_entrez);
    affy2entrez(i, inds_current_affy) = true;
end
    


save('affy_2_entrez.mat', 'affy2entrez','uniq_entrez','uniq_affy');
figure;hist(sum(affy2entrez,1),-1:200); xlabel('number of matches in affy'); ylabel('entrez count');
title(sprintf('%g%% of entrez have a unique match', sum( sum(affy2entrez,1) ==1) /size(affy2entrez,2) *100 ));
figure;hist(sum(affy2entrez,2),-1:200); xlabel('number of matches in entrez'); ylabel('affy count');
title(sprintf('%g%% of affy have a unique match', sum( sum(affy2entrez,2) ==1) /size(affy2entrez,1) *100 ));




z = textscan( fopen('affy_to_symbols.csv'), '%q %q','delimiter',',','HeaderLines',1);

affy = z{1};
symbol = z{2};
uniq_affy = unique(affy);
uniq_symbols = unique(symbol);
affy2symbol = false(length(uniq_affy), length(uniq_symbols));

for i=1:length(uniq_affy)
    currrent_affy = uniq_affy{i};
    inds_current_affy = strcmp(affy, currrent_affy);
    current_symbol = symbol(inds_current_affy);
    
    [~ ,inds_current_affy] = ismember(current_symbol, uniq_symbols);
    affy2symbol(i, inds_current_affy) = true;
end
    


save('affy_2_symbols.mat', 'affy2symbol','uniq_symbols','uniq_affy')



subplot(2,1,2);hist(sum(affy2symbol,1),100);
subplot(2,1,1);hist(sum(affy2entrez,1),100);
hist_affy2entrez = hist(sum(affy2entrez,1))
hist_affy2symbol = hist(sum(affy2symbol,1))