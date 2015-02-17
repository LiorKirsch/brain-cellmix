function [found_translation, entrez_or_symbol] = translateAffyToEntrez(list_of_affy_ids, entrez_or_symbol, get_only_1_1_match, entrez_or_symbol_subset)
% translate the list of affy ids to entrez
% if get_only_1_1_match is true it returns only translation that are 1-1
% if it is false it returns the first entrez id it finds.
% If you specify any values in "entrez_or_symbol_subset" this limits the
% returned values only to the values within this set.

switch entrez_or_symbol
    case 'symbol'
        tmp = load('affy_2_symbols.mat');
        target_id = tmp.uniq_symbols;
        afft2target = tmp.affy2symbol;
        uniq_affy = tmp.uniq_affy;
    case 'entrez'
        tmp = load('affy_2_entrez.mat');
        target_id = tmp.uniq_entrez;
        afft2target = tmp.affy2entrez;
        uniq_affy = tmp.uniq_affy;
end

if exist('entrez_or_symbol_subset','var')
    is_member = ismember(target_id, entrez_or_symbol_subset);
    target_id = target_id(is_member);
    afft2target = afft2target(:, is_member);
end

if get_only_1_1_match
    one2one_mask = sum(afft2target,2) == 1;
    uniq_affy = uniq_affy(one2one_mask);
    afft2target = afft2target(one2one_mask,:);
end

[found_translation, inds] = ismember(list_of_affy_ids, uniq_affy);


inds = inds(found_translation);
afft2target_subset = afft2target(inds,:);

if get_only_1_1_match
    assert( all(sum(afft2target_subset,2) == 1), 'should take only 1-1 matching genes');
end

[~,ind_in_entrez] = max(afft2target_subset,[],2);
entrez_or_symbol = target_id( ind_in_entrez );

    
end
