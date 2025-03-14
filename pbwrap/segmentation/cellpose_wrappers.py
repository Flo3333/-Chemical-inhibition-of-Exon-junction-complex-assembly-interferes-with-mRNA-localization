from .utils import unstack_slices, merge_channels, merge_channels_fromlists
from .utils import get_histogramm_highest_varation_value
from .utils import from2Dlabel_to3Dlabel

import numpy as np
import bigfish.stack as stack
import bigfish.segmentation as seg
import cellpose.models as models

from bigfish.stack import check_array, check_parameter
from pbwrap.integrity import check_sameshape
from ..errors import PbodySegmentationError, CellnumberError
from skimage.transform import resize


###### Image segmentation
def Nucleus_segmentation(dapi, diameter= 150, anisotropy= 3, use_gpu= False, use_3D_cellpose= False, model_type= 'Nuc_slices1.0', min_cell_number = 0) :
    """3D Nucleus segmentation using Cellpose from a dapi 3D grayscale image.

    Parameters
    ----------

        dapi :      np.ndarray, ndim = 3 (z,y,x). 
            The dapi should only contains the data to be analysed, prepropressing and out of focus filtering should be done prior to this operation. 
        diameter :  Int. 
            Average diameter of a nucleus in pixel. Used to rescale cellpose trained model to incoming data.
        anisotropy: Int. 
            Define the ratio between the plane (xy) resolution to the height (z) resolution. For a voxel size (300,100,100) use 3.
        use_gpu :   Bool. 
            Enable Cellpose build-in option to use GPU.
        min cell number : int.
            Raises CellNumberError(SegmentattionError) if less than min_cellnumber cells are segmented
                
    Returns
    -------
    
        Nucleus_label : np.ndarray 
            With same shape as dapi. Each object has a different pixel value and 0 is background.
    """

    #Integrity checks
    check_array(dapi, ndim= [2,3], dtype= [np.uint8, np.uint16, np.int32, np.int64, np.float32, np.float64])
    check_parameter(diameter= (int, float), anisotropy= (int), use_gpu= (bool))

    ndim = dapi.ndim

    #Segmentation
    nucleus_model = models.CellposeModel(gpu= use_gpu, model_type = model_type)
    channels = [0,0]

    if ndim == 3 :
        if use_3D_cellpose : 
            nucleus_label = nucleus_model.eval(dapi, diameter= diameter, channels = channels, anisotropy= anisotropy, do_3D= True)[0].astype(np.int64)
        else :
            list_dapi = unstack_slices(dapi)
            nucleus_label = nucleus_model.eval(list_dapi, diameter= diameter, channels= channels)[0]
            nucleus_label = np.array(nucleus_label, dtype = np.int64)
    
    else :
        nucleus_label = nucleus_model.eval(dapi, diameter= diameter, channels = channels)[0].astype(np.int64)
        nucleus_label = np.array(nucleus_label, dtype = np.int64)

    nucleus_label = seg.remove_disjoint(nucleus_label)

    if len(nucleus_label) < min_cell_number : raise CellnumberError("{0} nucleus were segmented, minimum cells number was set at {1}".format(len(nucleus_label), min_cell_number))
    return nucleus_label




def Cytoplasm_segmentation(cy3, dapi= None, diameter= 250, maximal_distance= 100, use_gpu= False, model_type = "Hek_2.0", min_cell_number= 0) :
    
    """Due to bad performance using 3D cellpose with cy3 channel image. A 2D cell segmentation is performed for each slice and a 3D labelling is performed using a closest centroid method.

    Parameters
    ----------

        cy3 :      np.ndarray, ndim = 3 (z,y,x). 
            cy3 should only contains the data to be analysed, prepropressing and out of focus filtering should be done prior to this operation. 
        diameter :  Int. 
            Average diameter of a cell in pixel. Used to rescale cellpose trained model to incoming data.
        use_gpu :   Bool. 
            Enable Cellpose build-in option to use GPU.
    
    Returns
    -------
    
        cytoplasm_labels : List[np.ndarray] 
            List of numpy arrays with shape(x,y). Eache element correspond to a z plane and each object has a different pixel value and 0 is background.
    """
    #Integrity checks
    check_array(cy3, ndim= [2,3], dtype= [np.uint8, np.uint16, np.int32, np.int64, np.float32, np.float64])
    check_parameter(diameter= (int, float), dapi= (np.ndarray, type(None)), use_gpu= (bool))

    #Segmentation
    if cy3.ndim == 3 : cytoplasm_slices = unstack_slices(cy3)
    cytoplasm_model = models.CellposeModel(gpu= use_gpu, model_type = model_type)
 

    if type(dapi) == type(None) : 
        channels = [0,0]
        if cy3.ndim == 3 : image = cytoplasm_slices 
        else : image = cy3
    else :
        check_sameshape(cy3, dapi)
        channels = [1,2]
        if dapi.ndim == 3 :
            nucleus_slices = unstack_slices(dapi)
            image = merge_channels_fromlists(cytoplasm_slices, nucleus_slices)
        else :
            image = merge_channels(cy3,dapi)


    cytoplasm_labels = cytoplasm_model.eval(image, diameter= diameter, channels= channels, do_3D= False)[0]
   
    return cytoplasm_labels

def Cytoplasm_segmentation_old(cy3, dapi= None, diameter= 250, maximal_distance= 100, use_gpu= False) :

    """Due to bad performance using 3D cellpose with cy3 channel image. A 2D cell segmentation is performed for each slice and a 3D labelling is performed using a closest centroid method.

    Parameters
    ----------

        cy3 :      np.ndarray, ndim = 3 (z,y,x). 
            cy3 should only contains the data to be analysed, prepropressing and out of focus filtering should be done prior to this operation. 
        diameter :  Int. 
            Average diameter of a cell in pixel. Used to rescale cellpose trained model to incoming data.
        use_gpu :   Bool. 
            Enable Cellpose build-in option to use GPU.
    
    Returns
    -------
    
        cytoplasm_labels : List[np.ndarray] 
            List of numpy arrays with shape(x,y). Eache element correspond to a z plane and each object has a different pixel value and 0 is background.
    """
    #Integrity checks
    check_array(cy3, ndim= [2,3], dtype= [np.uint8, np.uint16, np.int32, np.int64, np.float32, np.float64])
    check_parameter(diameter= (int), dapi= (np.ndarray, type(None)), use_gpu= (bool))


    #image downscale
    scale_factor = diameter / 30
    if cy3.ndim == 2 :
        cy3_rescaled = resize(cy3, (round(cy3.shape[0] / scale_factor), round(cy3.shape[1] / scale_factor)), anti_aliasing= True)
        
    else : 
        cy3_rescaled = resize(cy3, (round(cy3.shape[0] / scale_factor), round(cy3.shape[1] / scale_factor), round(cy3.shape[2] / scale_factor)), anti_aliasing= True)
        

    #Segmentation
    min_objct_size = int(round((np.pi * (diameter/2)**2) /5)) # area in pixel
    if cy3.ndim == 3 : cytoplasm_slices = unstack_slices(cy3_rescaled)
    cytoplasm_model = models.Cellpose(gpu= use_gpu, model_type = "cyto")
 

    if type(dapi) == type(None) : 
        channels = [0,0]
        if cy3.ndim == 3 : image = cytoplasm_slices 
        else : image = cy3_rescaled
    else :
        check_sameshape(cy3, dapi)
        channels = [1,2]
        if dapi.ndim == 3 :
            dapi_rescaled = resize(dapi, (round(dapi.shape[0] / scale_factor), round(dapi.shape[1] / scale_factor), round(dapi.shape[2] / scale_factor)), anti_aliasing= True)
            nucleus_slices = unstack_slices(dapi_rescaled)
            image = merge_channels_fromlists(cytoplasm_slices, nucleus_slices)
        else :
            dapi_rescaled = resize(dapi, (round(dapi.shape[0] / scale_factor), round(dapi.shape[1] / scale_factor)), anti_aliasing= True)
            image = merge_channels(cy3_rescaled,dapi_rescaled)


    cytoplasm_labels = cytoplasm_model.eval(image, diameter= 30, channels= channels, do_3D= False)[0]
    cytoplasm_labels = resize(cytoplasm_labels, cy3.shape, preserve_range= True)
    cytoplasm_labels = np.array(cytoplasm_labels, dtype= np.int64)


    if cy3.ndim == 3 : 
        for z in range(0,len(cytoplasm_labels)): cytoplasm_labels[z] = seg.clean_segmentation(cytoplasm_labels[z], small_object_size= min_objct_size, delimit_instance=True, smoothness= 3,  fill_holes= True)
        cytoplasm_label = from2Dlabel_to3Dlabel(cytoplasm_labels, maximal_distance= maximal_distance)
    else : 
        cytoplasm_label = seg.clean_segmentation(cytoplasm_labels, small_object_size= min_objct_size, delimit_instance=True, fill_holes= True)

    cytoplasm_label = seg.remove_disjoint(cytoplasm_label)
    cytoplasm_label[cytoplasm_label < 1] = 0
    cytoplasm_label[cytoplasm_label > 0] = 1
    cytoplasm_label= seg.label_instances(np.array(cytoplasm_label, dtype = bool)) 

    return cytoplasm_label





def pbody_segmentation(egfp, sigma = 2, threshold= None, thresh_penalty= 1, small_object_size= None, fill_holes= True) :
    """Performs Pbody segmentation on 2D or 3D egfp numpy array.
        Apply a log filter to the image which kernel is defined by sigma.
        Then a threshold is applied, if none is given compute automatic threshold from highest variation point in array histogram.
        if peaks_min_distance other than 0 is given, performs a watershed segmentation.
    
    Parameters
    ----------
        egfp : np.ndarray(y,x) or (z,y,x)
        beta : scalar
            Parameter representing the difficulty to cross region with high intensity gradient.
            
    Returns
    -------
        Pbody_label : np.ndarray(y,x) or (z,y,x)
    """

    check_parameter(egfp = (np.ndarray))
    check_array(egfp, ndim= [2,3])
    dim = egfp.ndim
    #Segmentation
    mask = stack.log_filter(egfp,sigma)
    if dim == 2 :
        if threshold == None : 
            threshold = get_histogramm_highest_varation_value(mask)
            threshold *= thresh_penalty
        mask = seg.thresholding(mask, threshold)
        mask = seg.clean_segmentation(mask, small_object_size=small_object_size, fill_holes=fill_holes)
    else : 
        for z in range(0,egfp.shape[0]):
            if threshold == None : 
                threshold = get_histogramm_highest_varation_value(mask[z])
                threshold *= thresh_penalty
            mask[z] = seg.thresholding(mask[z], threshold)
            mask[z] = seg.clean_segmentation(np.array(mask[z],dtype= bool), small_object_size=small_object_size, fill_holes=fill_holes)

    #labelling
    if dim == 2 : 
        egfp_label = seg.label_instances(mask)
        egfp_label = np.array(egfp_label, dtype= np.int64)

    else :
        slice_list = unstack_slices(mask)
        egfp_label = from2Dlabel_to3Dlabel(slice_list, maximal_distance= 10)



    if len(egfp_label) == 0 : raise PbodySegmentationError("No pbody was segmentated.")
    return egfp_label