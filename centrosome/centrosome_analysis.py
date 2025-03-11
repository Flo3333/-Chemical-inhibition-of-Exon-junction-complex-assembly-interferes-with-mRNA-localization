import os
import CustomPandasFramework.centrosome.analysis as analysis
import CustomPandasFramework.centrosome.cellular_cycle as cycle_analysis
from bigfish.classification import get_features_name


##############################################
##############   PARAMETERS     ##############
input_path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/output/20240418_15-14-14/result_tables/'
rna_list = None
# rna_list = ['AKAP9', 'C-','DYNC1H1', 'FLAP XY', 'GAS8', 'GOLGA4', 'IQCB1', 'KIAA0753', 'KIF1C', 'NIN', 'NUMA1', 'PCNT', 'SSH1'] #stack_EJC_HeLacentrin #CHECK2 deleted
rna_list = ['AKAP9', 'C-','FLAP XY', 'KIF1C', 'NIN', 'NUMA1', 'PCNT',] #stack_EJC_Hek4Y
index = ['rna', 'treatment', 'well', 'fov']

frameon = True

do_Excel_extracts = True
extract_columns = ['proportion_rna_in_nuc','nb_rna_out_nuc','nb_rna_in_nuc','proportion_nuc_area','cell_area','nuc_area','centrosome_number','index_mean_distance_centrosome','index_median_distance_centrosome','index_rna_centrosome','proportion_rna_centrosome','index_centrosome_dispersion']

do_individual_box_plots = True
box_plots_columns = ['proportion_rna_in_nuc','nb_rna_out_nuc','nb_rna_in_nuc','index_mean_distance_centrosome','index_median_distance_centrosome','proportion_rna_centrosome']

do_rna_number_per_cell_plot = True

do_rna_distribution_plot = True
rna_distribution_list = ['proportion_rna_centrosome'] #if None does all

do_overall_bar_plots = True
overal_bar_plots_list = ['proportion_rna_centrosome', 'index_median_distance_centrosome']

do_outer_nuclei_caracterisation = False
outer_nuclei_rna_list = ['AKAP9', 'GOLGA4', 'DYNC1H1', 'NIN', 'KIF1C']
# outer_nuclei_measures = ["index_median_distance_cell"]
outer_nuclei_measures = get_features_name(names_features_distance= True)

do_focci_analysis = False
focci_analysis_measures = ['cluster_number', 'nucleus_cluster_number', 'clustered_spots_number', 'unclustered_spots_number', 'clustered_spots_fraction']
focci_rna_list = ['DYNC1H1']

do_Tukey_hsd_tiles= True
Tukey_rna_list = rna_list
Tukey_measures_list = ['clustered_spots_fraction', 'cluster_number', 'clustered_spots_number', 'index_median_distance_centrosome', 'proportion_rna_centrosome']

statistical_test_significance = 0.01

do_cellular_cycle_analysis = False

################################################
################################################
#                   SCRIPT                     #
print("Centrosome analysis starts.\n")

#Prepare result dir
if not input_path.endswith('/') : input_path += '/'
output_path = input_path.replace('result_tables', 'result_plots')
os.makedirs(output_path, exist_ok= True)

#Analysis

if do_Excel_extracts : 
    print("Extracting Excel results...\n")
    analysis.Cell_extract(input_path, output_path, extract_columns=extract_columns, index_keys=index + ['id'], rna_list= rna_list)

if do_individual_box_plots : 
    print("Starting individual box plots :")
    analysis.individual_box_plots(input_path, output_path, box_plots_columns, rna_list=rna_list, index_keys= index + ['id'])
    print("")

if do_rna_number_per_cell_plot :
    print("Rna per cell plot...\n")
    analysis.individual_violin_rna_number_per_cell(input_path, output_path, rna_list=rna_list, frameon=frameon)
    analysis.rna_number_per_cell_plot(input_path, output_path, rna_list=rna_list, frameon=frameon)

if do_rna_distribution_plot :
    print("Starting rna distribution plots :")
    analysis.cell_distribution(input_path, output_path, measure_list=rna_distribution_list, rna_list=rna_list,frameon=frameon )

if do_overall_bar_plots :
    print("Starting overall plots :")
    analysis.overall_plot(input_path, output_path, measures= overal_bar_plots_list, rna_list=rna_list, frameon=frameon)
    analysis.p_value(input_path, output_path, measure_list= overal_bar_plots_list, significance=statistical_test_significance, rna_list=rna_list, frameon=frameon, folder_name= 'overall_plots')

if do_outer_nuclei_caracterisation :
    print("Starting outer nuclei caracterisation :")
    analysis.overall_plot(input_path, output_path, measures= outer_nuclei_measures, rna_list= outer_nuclei_rna_list, folder_name= "outer_nuclei_caracterisation", frameon=frameon)
    analysis.p_value(input_path, output_path, measure_list= outer_nuclei_measures, significance=statistical_test_significance, rna_list=outer_nuclei_rna_list, frameon=frameon, folder_name="outer_nuclei_caracterisation")

if do_focci_analysis :
    print("Starting focci analysis : ")
    analysis.individual_violin_plots(input_path, output_path, measures= focci_analysis_measures, rna_list= focci_rna_list, folder_name= "focci_analysis", frameon=frameon)
    analysis.p_value(input_path, output_path, measure_list= focci_analysis_measures, significance=statistical_test_significance, rna_list=focci_rna_list, frameon=frameon, folder_name="focci_analysis")

if do_Tukey_hsd_tiles :
    print("Starting Tukey-hsd tiles :")
    analysis.Tukey_hsd_tiles(input_path, output_path, measure_list=Tukey_measures_list, rna_list=Tukey_rna_list)

if do_cellular_cycle_analysis :
    print("Starting cellular cycle analysis :")