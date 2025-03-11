import os, warnings
import numpy as np
import CustomPandasFramework.centrosome as framework
import CustomPandasFramework.log as log
import CustomPandasFramework.computer_interface as interface
import pbwrap.detection as detection
import pbwrap.segmentation as segmentation
import pbwrap.quantification as quant
import pbwrap.plot as plot
import bigfish.stack as stack
import bigfish.plot as seg_plot
from CustomPandasFramework.operations import add_data
from bigfish.multistack import match_nuc_cell
from bigfish.detection import detect_spots
from numpy import NaN

path_in = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/input'
path_out = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/output'

#Full list stack_EJC_HeLacentrin : ['AKAP9', 'C-', 'CHEK2', 'DYNC1H1', 'FLAP XY', 'GAS8', 'GOLGA4', 'IQCB1', 'KIAA0753', 'KIF1C', 'NIN', 'NUMA1', 'PCNT', 'SSH1']
# 1.2 penalty : ['AKAP9', 'C-', 'CHEK2', 'FLAP XY', 'KIAA0753', 'KIF1C', 'NIN', 'NUMA1', 'PCNT']
# 0.6 penalty : ['GAS8', 'GOLGA4', 'IQCB1', 'SSH1', 'DYNC1H1'] --> DYNC1H1 and GOLGA4 to lower even more
# 0.4 penalty : ['DYNC1H1']

rna_penalty_dict = {}
rna_penalty_dict.update(dict.fromkeys([''], 0.1))
rna_penalty_dict.update(dict.fromkeys(['NIN','GAS8', 'GOLGA4', 'SSH1','PCNT','DYNC1H1'], 0.3))
rna_penalty_dict.update(dict.fromkeys(['IQCB1'], 0.4))
rna_penalty_dict.update(dict.fromkeys(['AKAP9',], 0.6))
rna_penalty_dict.update(dict.fromkeys(['C-', 'CHEK2', 'FLAP XY', 'KIAA0753', 'KIF1C', 'NUMA1'], 1.2))

remove_artifact_dict = dict.fromkeys(['DYNC1H1'], False)


## Parameters
nucleus_channel = 'DAPI'
rna_channel = 'Cy3'
centrosome_channel = 'Alexa 647'
gene_list = None#
voxel_size = (300, 103, 103)
spot_multiplicator = 0.75
spot_size = (350*spot_multiplicator, 150*spot_multiplicator, 150*spot_multiplicator)
crop_z = (2,15)

#Detection
RNA_THRESHOLD = None
z_stacks_threshold_computation = (3,7)
rna_penalty_default = 1
show_progression = True

#Deconvolution
rna_alpha = 0.5
rna_beta = 1
deconvolution_timer = 60 #s
cluster_size = 600
cluster_min_spot_number = 3

# Segmentation
nuc_model = "nuclei"
cell_model = "cyto2"
nucleus_size = 70 #px - ~diameter
cell_size = 100 #px - ~ diameter
centrosome_size = (400,300,300) #nm - diameter
centrosome_penalty_threshold = 2
##Choose one
Log_method = False
Dog_method = True

#Notification
enable_notifications = False


##################################################  Pipeline ###############################################################
print("Starting analysis for genes : ", gene_list)
print("Threshold is set to {0}".format(RNA_THRESHOLD))
print("Default threshold penalty is set to {0}".format(rna_penalty_default))

path_out = interface.add_datetime_folder(path_out)
channel = [nucleus_channel, rna_channel, centrosome_channel]

os.makedirs(path_out)
run_log = log.run_log('log.txt', path=path_out)
rna_computed = []
parameters_log = log.parameter_log('parameters.txt', path=path_out)
parameters_log.add_parameters(channel, gene_list, voxel_size, spot_size, RNA_THRESHOLD, rna_penalty_default, rna_penalty_dict, show_progression, rna_alpha, rna_beta, deconvolution_timer, nuc_model, cell_model, nucleus_size, cell_size)
parameters_log.write()

if not path_out.endswith('/') : path_out += '/'
tables_path = path_out + 'result_tables/'
os.makedirs(tables_path, exist_ok= True)

Input = framework.make_Input(path_in=path_in, channel=channel)
Acquisition = framework.empty_frame()
Cell = framework.empty_frame()
if type(gene_list) != type(None) : 
    Input_idx = Input.query('rna_name in {0}'.format(gene_list)).index
    Input = Input.loc[Input_idx,:]

#analysis loop
for rna in Input.index.get_level_values(0).unique() :
    print("\ncomputing rna : {0}".format(rna))
    #rna setting
    rna_penalty = rna_penalty_dict.setdefault(rna, rna_penalty_default)
    do_remove_artifact = remove_artifact_dict.setdefault(rna, True)
    print('remove artifact : {0}'.format(do_remove_artifact))

    #Spot detection
    rna_list = framework.open_image(Input, path_in, channel= rna_channel, rna_name= rna, output_list= True, crop_z=crop_z)
    centrosome_list = framework.open_image(Input, path_in, channel= centrosome_channel, rna_name= rna, output_list= True, crop_z=crop_z)
    print('image number : ', len(rna_list))
    print("computing threshold...")
    print("penalty is set to {0}".format(rna_penalty))
    rna_threshold = RNA_THRESHOLD
    if type(rna_threshold) == type(None) : rna_threshold = detection.compute_auto_threshold(rna_list, voxel_size=voxel_size, spot_radius=spot_size, im_number= 10, crop_zstack=z_stacks_threshold_computation) * rna_penalty
    centrosome_threshold = detection.compute_auto_threshold(centrosome_list, voxel_size=voxel_size, spot_radius=spot_size, im_number= 10, crop_zstack=z_stacks_threshold_computation) * centrosome_penalty_threshold
    print('threshold set to : ', rna_threshold)

    dapi_gen = framework.open_image(Input, path_in, nucleus_channel, rna_name= rna, crop_z= crop_z)
    cy3_gen = framework.open_image(Input, path_in, rna_channel, rna_name= rna, crop_z=crop_z)
    egfp_gen = framework.open_image(Input, path_in, centrosome_channel, rna_name= rna, crop_z=crop_z)

    #fov loop
    rna_Input = Input.loc[rna,:]
    print("RNA_Input\n", rna_Input)
    for idx in rna_Input.reset_index(level=3).index.unique() :
        treatment = idx[0] if idx[0] != NaN else 'untreated'
        well = idx[1]
        fov = idx[2]
        cell_line = rna_Input.at[(treatment, well, fov, centrosome_channel), 'cell_line']
        fov_id = rna_Input.at[(treatment, well, fov, centrosome_channel),'id']
        print( "\ncomputing fov {0}.".format(idx))
        run_log.update('fov {0}'.format(idx), rna_computed)

        #Pre-processing
        print("  preprocessing fov...")
        dapi = next(dapi_gen)
        cy3 = next(cy3_gen)
        egfp = next(egfp_gen)
        cy3_max_proj = stack.maximum_projection(cy3)
        cy3_mean_proj = stack.mean_projection(cy3)
        dapi_mean_proj = stack.mean_projection(dapi)
        egfp_max_proj = stack.maximum_projection(egfp)
        egfp_mean_proj = stack.mean_projection(egfp)

        #Spot detection
        print("  spot detection...")
        rna_spots = detect_spots(cy3, threshold= rna_threshold, voxel_size=voxel_size, spot_radius=spot_size)
        centrosome_spots = detect_spots(egfp, threshold= centrosome_threshold, voxel_size=voxel_size, spot_radius=spot_size)

        #segmentation
        print("  running segmentations...")
        nucleus_label = segmentation.Nucleus_segmentation(dapi_mean_proj, diameter= nucleus_size, use_gpu= True, model_type= nuc_model)
        cell_label = segmentation.Cytoplasm_segmentation(cy3_mean_proj, dapi_mean_proj, diameter= cell_size, use_gpu= True, model_type= cell_model)
        # centrosome_label = segmentation.centrosome_segmentation_candidate_regions(egfp, centrosome_size=centrosome_size, voxel_size=voxel_size, gaussian_fit= 'FWHM', threshold_penalty=centrosome_penalty_threshold, DoG_filter_method=Dog_method, LoG_filter_method=Log_method)
        nucleus_label, cell_label = match_nuc_cell(nuc_label=nucleus_label, cell_label= cell_label, single_nuc=True, cell_alone=False)

        #deconvolution / clustering
        print("  clusters deconvolution and scanning...")
        deconvoluted_rna_spots = detection.cluster_deconvolution(cy3, rna_spots, spot_radius= spot_size, voxel_size= voxel_size, alpha= rna_alpha, beta= rna_beta, timer=deconvolution_timer)
        if do_remove_artifact : deconvoluted_rna_spots = detection.remove_artifact(deconvoluted_rna_spots, 103*10, voxel_size=voxel_size, spot_density= 2)
        clustering_dict = detection.cluster_detection(deconvoluted_rna_spots, voxel_size= voxel_size, radius= cluster_size, nb_min_spots= cluster_min_spot_number, keys_to_compute=['clusters_dataframe', 'clustered_spots_dataframe'])
        Clusters_fov = clustering_dict['clusters_dataframe']
        Spots_fov = clustering_dict['clustered_spots_dataframe']
        cluster_coords = detection.get_centroids_array(Clusters_fov)
        clustered_spots = detection.get_centroids_array(Spots_fov.query("not cluster_id.isna()"))
        unclustered_spots = detection.get_centroids_array(Spots_fov.query("cluster_id.isna()"))
        del clustering_dict
        if len(deconvoluted_rna_spots) == 0 : 
            print("\033[91mNo spot were detected :\033[0m computing next fov.")
            run_log.failure()
            continue

        #Controls visuals
        print("  saving visuals...")
        segmentation_path = path_out + "control_visuals/segmentation/{0}/{1}/".format(rna, treatment)
        os.makedirs(segmentation_path, exist_ok=True)
        seg_plot.plot_segmentation_boundary(cy3_mean_proj, cell_label= cell_label, nuc_label= nucleus_label, boundary_size= 2, remove_frame=True, contrast= True, path_output= segmentation_path + "cell_segmentation_{0}".format(fov), show=False)
        detection_path = path_out + "control_visuals/detection/{0}/{1}/".format(rna, treatment)
        os.makedirs(detection_path, exist_ok=True)
        plot.output_spot_tiffvisual(cy3_max_proj, [rna_spots, deconvoluted_rna_spots, clustered_spots, unclustered_spots], detection_path + "rna_spot_detection_{1}_{0}".format(fov, well), dot_size=2)

        #Quantification
        print("  starting fov quantification ...")
        other_images = {"smfish" : cy3_max_proj, "centrosome_mip" : egfp_mean_proj}
        other_coords = {"clusters_coords" : cluster_coords, "clustered_spots" : clustered_spots, "unclustered_spots" : unclustered_spots, "centrosome_spots" : centrosome_spots}
        fov_results = quant.extract_cell(cell_label=cell_label, ndim= 3, image= cy3_max_proj,others_image=other_images, nuc_label=nucleus_label, rna_coord= deconvoluted_rna_spots, others_coord=other_coords)
        
        #Cell_loop
        count_cell_discarded = 0
        count_cell_computed = 0
        cell_path = path_out + "control_visuals/cells/{0}/{1}/{2}/".format(rna, treatment, fov_id)
        os.makedirs(cell_path, exist_ok=True)
        Cell_result = framework.empty_frame()

        warnings.simplefilter("ignore") 
        for cell_id, cell in enumerate(fov_results) :            
            try :
                Cell_result = framework.add_cell(Cell_df= Cell_result, acquisition_id= fov_id, cell= cell, dapi_stack= dapi, voxel_size=voxel_size)
            except quant.QuantificationError as error :
                if str(error) == 'No centrosome found' : count_cell_discarded += 1
                else: raise error
            else :
                framework.centrosome_cell_plot(cell, Cell_result.tail(1), path_output= cell_path)
                count_cell_computed += 1
            finally : 
                warnings.simplefilter("default")
        
        #centrosome final visual
        if "centrosome_coords" in Cell_result.columns :
            final_centrosome_spots = np.array(sum(list(Cell_result['centrosome_coords'].apply(list)), []), dtype= int)
            plot.output_spot_tiffvisual(egfp_max_proj, [centrosome_spots, final_centrosome_spots], detection_path + "centrosome_spot_detection_{1}_{0}".format(fov, well), dot_size=4)
            
        Cell = add_data(Cell, Cell_result)
        Acquisition = framework.add_fov(Acquisition, acquisition_id= fov_id, cell_line=cell_line, rna=rna, well= well, fov= fov, treatement=treatment, rna_threshold= rna_threshold, cell_computed=count_cell_computed, cell_discarded= count_cell_discarded)
        run_log.sucess()

    rna_computed += [str(rna)]
    Input.reset_index(drop=True).to_feather(tables_path + 'Input.feather')
    Acquisition.reset_index(drop=True).to_feather(tables_path + 'Acquisition.feather')
    Cell.reset_index(drop=True).to_feather(tables_path + 'Cell.feather')

print("pipeline end. Saving results...")

Input.reset_index(drop=True).to_excel(tables_path + 'Input.xlsx')
Acquisition.reset_index(drop=True).to_excel(tables_path + 'Acquisition.xlsx')
Cell.reset_index(drop=True).to_excel(tables_path + 'Cell.xlsx')

print("Done.")