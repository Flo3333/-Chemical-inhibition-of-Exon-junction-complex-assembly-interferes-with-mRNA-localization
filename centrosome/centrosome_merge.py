import CustomPandasFramework.operations.merge as merge

path_list = [ # (paths, rna_list)
    ('/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/output/set_stack_EJC_HeLacentrin_final_20240118/raw_results/20240111_11-06-58/result_tables', ['DYNC1H1']),
    ('/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/output/set_stack_EJC_HeLacentrin_final_20240118/raw_results/20240117_15-46-33/result_tables', ['AKAP9', 'C-', 'CHEK2', 'FLAP XY', 'GAS8', 'GOLGA4', 'IQCB1', 'KIAA0753', 'KIF1C', 'NIN', 'NUMA1', 'PCNT', 'SSH1'])
    ] 
save_path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/output/set_stack_EJC_HeLacentrin_final_20240118/final_result_tables/'
table_names = ['Acquisition', 'Cell']

tables = {}

for result_path, rna_list in path_list :
    new_tables = merge.open_tables(result_path, table_names, rna_list=rna_list)
    if len(tables) > 0 : duplicates = merge.look_for_duplicates(tables['Acquisition'], new_tables['Acquisition'], ['rna'])
    else : duplicates = []
    if len(duplicates) > 0 : raise Exception("Duplicate rna were found at {0}\n{1}".format(result_path, duplicates))

    tables = merge.merge_tables(tables, new_tables)

merge.save_merged_tables(tables, path_out=save_path)
