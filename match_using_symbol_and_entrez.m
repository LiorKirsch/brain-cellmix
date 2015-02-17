function match_using_symbol_and_entrez(symbols_A, entrez_A, symbols_B, entrez_B)

unique_entrez = unique( [entrez_A; entrez_B]);
unique_symbols = unique( [symbols_A; symbols_B]);

[~, indA_entrez] = ismember(entrez_A, unique_entrez);
[~, indB_entrez] = ismember(entrez_B, unique_entrez);
ismember_entrez_matrix_A = sparse(1:length(indA_entrez), indA_entrez, ones(length(indA_entrez),1),...
                            length(indA_entrez), length(unique_entrez) );
ismember_entrez_matrix_B = sparse(1:length(indB_entrez), indB_entrez, ones(length(indB_entrez),1),...
                            length(indB_entrez), length(unique_entrez) );
                        
                        
[~, indA_symbols] = ismember(entrez_A, unique_symbols);
[~, indB_symbols] = ismember(entrez_B, unique_symbols);
ismember_symbols_matrix_A = sparse(1:length(indA_symbols), indA_symbols, ones(length(indA_symbols),1),...
                            length(indA_symbols), length(unique_symbols) );
ismember_symbols_matrix_B = sparse(1:length(indB_symbols), indB_symbols, ones(length(indB_symbols),1),...
                            length(indB_symbols), length(unique_symbols) );

                        
end