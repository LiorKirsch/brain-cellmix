% findGenesInDifferentAnimal()


mouse_gene_entrez = {'20257', 	'108013', '319586', '76183'}';

[gene_to_group_mouse, gene_to_group_monkey, homologous_group_id, gene_id_A, gene_id_B] =  gene_to_homolog_group('mouse_laboratory','rhesus_macaque', mouse_gene_entrez, 'entrez_gene_ID');



mouse_gene_names = {'Stmn2', 'Celf4',  'Syt1', 'Gfap', 'Aqp4', 'Fgfr3', 'Slc1a2', 'Gjb6', 'Mbp', 'Sox10', 'Mag','Mog'}';

[gene_to_group_mouse, gene_to_group_monkey, homologous_group_id, gene_id_A, gene_id_B] =  gene_to_homolog_group('mouse_laboratory','rhesus_macaque', mouse_gene_names, 'symbol',[],'entrez_gene_ID');



{'';;'SOX10';}

MBP - LOC720908
SYT1 - Syt1
Mag - MAG
STMN2 - STMN2
GFAP - GFAP
MOG - MOG
FGFR3 - FGFR3
GJB6 - GJB6
SLC1A2 - SLC1A2