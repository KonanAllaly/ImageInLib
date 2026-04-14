#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h> 
#include "common_functions.h"
#include "../src/vtk_params.h"
#include "common_vtk.h"
#include "nifti1_io.h"
    
    typedef struct
    {
        //Dimensions
        size_t length, width, heigth;
        
        //Spacing
        dataType sx, sy, sz;
        
        //Origin
        float ox, oy, oz;
        
        //Data
        int datatype;          /* NIfTI datatype code (e.g. DT_INT16)   */
        int bytes_per_voxel;   /* e.g. 2 for INT16, 4 for FLOAT32       */
        void** data;           /* data[z] points to slice z             */

    } NiiContainer;
    
    bool free_Container(NiiContainer* imageDataStr);

    dataType get_voxel(const NiiContainer* imageDatastr, size_t x, size_t y, size_t z);
    
    bool load_nii(const char* filepath, NiiContainer* vol);
    
    bool saveNiiVolumeAsVtk(const NiiContainer* imageDatastr, const char* vtk_path, vtkDataForm dataForm);

#ifdef __cplusplus
}
#endif