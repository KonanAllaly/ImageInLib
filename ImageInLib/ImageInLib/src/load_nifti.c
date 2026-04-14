#include "load_nifti.h"

bool free_Container(NiiContainer* imageDatastr)
{
    if (!imageDatastr)
    {
        fprintf(stderr, "[free_Container] ERROR: NULL pointer provided\n");
		return false;
    }

    if (imageDatastr->data)
    {
        //Free each slice, then the pointer array
        for (size_t k = 0; k < imageDatastr->heigth; k++)
        {
            free(imageDatastr->data[k]);
        }
        free(imageDatastr->data);
        imageDatastr->data = NULL;
    }
	return true;
}

dataType get_voxel(const NiiContainer* imageDatastr, size_t x, size_t y, size_t z)
{
    size_t xd = x_new(x, y, imageDatastr->length);

    switch (imageDatastr->datatype)
    {
    case DT_UINT8:
        return (dataType)((uint8_t*)imageDatastr->data[z])[xd];
    case DT_INT8:
        return (dataType)((int8_t*)imageDatastr->data[z])[xd];
    case DT_UINT16:
        return (dataType)((uint16_t*)imageDatastr->data[z])[xd];
    case DT_INT16:
        return (dataType)((int16_t*)imageDatastr->data[z])[xd];
    case DT_UINT32:
        return (dataType)((uint32_t*)imageDatastr->data[z])[xd];
    case DT_INT32:
        return (dataType)((int32_t*)imageDatastr->data[z])[xd];
    case DT_FLOAT32:
        return ((dataType*)imageDatastr->data[z])[xd];
    case DT_FLOAT64:
        return (dataType)((double*)imageDatastr->data[z])[xd];
    default:
        fprintf(stderr, "get_voxel: unsupported datatype %d\n", imageDatastr->datatype);
        return 0.0f;
    }
}

bool load_nii(const char* filepath, NiiContainer* imageDatastr)
{
    if (!filepath || !imageDatastr) 
    {
		fprintf(stderr, "[load_nii] ERROR: NULL pointer provided for filepath or imageDatastr\n");
        return false;
    } 

    // Read header and data
    nifti_image* nim = nifti_image_read(filepath, 1);
    if (!nim)
    {
        fprintf(stderr, "[load_nii] ERROR: cannot open '%s'\n", filepath);
        return false;
    }

    if (nim->ndim < 3)
    {
        fprintf(stderr, "[load_nii] ERROR: image has only %d dimension(s)\n",
            nim->ndim);
        nifti_image_free(nim);
        return false;
    }

    //Dimensions
    imageDatastr->length = nim->nx;
    imageDatastr->width = nim->ny;
    imageDatastr->heigth = nim->nz;

    //Spacing
    imageDatastr->sx = nim->dx;
    imageDatastr->sy = nim->dy;
    imageDatastr->sz = nim->dz;

    //Origin
    if (nim->sform_code > 0)
    {
        imageDatastr->ox = nim->sto_xyz.m[0][3];
        imageDatastr->oy = nim->sto_xyz.m[1][3];
        imageDatastr->oz = nim->sto_xyz.m[2][3];
    }
    else if (nim->qform_code > 0)
    {
        imageDatastr->ox = nim->qto_xyz.m[0][3];
        imageDatastr->oy = nim->qto_xyz.m[1][3];
        imageDatastr->oz = nim->qto_xyz.m[2][3];
    }
    else
    {
        //No transform available — origin unknown
        imageDatastr->ox = imageDatastr->oy = imageDatastr->oz = 0.0f;
    }

    //Datatype
    imageDatastr->datatype = nim->datatype;
    imageDatastr->bytes_per_voxel = nim->nbyper;

    //Allocate data[nz] — one pointer per slice
    size_t slice_bytes = (size_t)imageDatastr->length * imageDatastr->width * imageDatastr->bytes_per_voxel;

    imageDatastr->data = (void**)malloc(imageDatastr->heigth * sizeof(void*));
    if (!imageDatastr->data)
    {
        fprintf(stderr, "[load_nii] ERROR: malloc failed for slice pointers\n");
        nifti_image_free(nim);
        return false;
    }

    for (size_t k = 0; k < imageDatastr->heigth; k++)
    {
        imageDatastr->data[k] = malloc(slice_bytes);
        if (!imageDatastr->data[k])
        {
            fprintf(stderr, "[load_nii] ERROR: malloc failed for slice %d\n", k);
            
            //Free already allocated slices
            for (size_t i = 0; i < k; i++) free(imageDatastr->data[i]);
            free(imageDatastr->data);
            imageDatastr->data = NULL;
            nifti_image_free(nim);
            return false;
        }
        //Copy slice k from the flat nifti buffer into data[k]
        memcpy(imageDatastr->data[k], (char*)nim->data + k * slice_bytes, slice_bytes);
    }

    nifti_image_free(nim);
    return true;
}

bool saveNiiVolumeAsVtk(const NiiContainer* imageDatastr, const char* vtk_path, vtkDataForm dataForm)
{
    if (!imageDatastr || !vtk_path)
    {
        return false;
    } 
    
    const size_t points_in_slice = (size_t)imageDatastr->length * imageDatastr->width;

    //Fill Vtk_File_Info metadata ──────────────────────
    Vtk_File_Info vtkMetaInfo;
    memset(&vtkMetaInfo, 0, sizeof(Vtk_File_Info));

    vtkMetaInfo.dimensions[0] = imageDatastr->length;
    vtkMetaInfo.dimensions[1] = imageDatastr->width;
    vtkMetaInfo.dimensions[2] = imageDatastr->heigth;

    vtkMetaInfo.spacing[0] = (dataType)imageDatastr->sx;
    vtkMetaInfo.spacing[1] = (dataType)imageDatastr->sy;
    vtkMetaInfo.spacing[2] = (dataType)imageDatastr->sz;

    vtkMetaInfo.origin[0] = (dataType)imageDatastr->ox;
    vtkMetaInfo.origin[1] = (dataType)imageDatastr->oy;
    vtkMetaInfo.origin[2] = (dataType)imageDatastr->oz;

    vtkMetaInfo.vDataType = dta_Flt;
    vtkMetaInfo.operation = copyFrom;   //fillPtr will read FROM dataPointer

    vtkMetaInfo.dataPointer = (dataType**)malloc(imageDatastr->heigth * sizeof(dataType*));
    if (!vtkMetaInfo.dataPointer) 
     {
         fprintf(stderr, "[saveNiiVolumeAsVtk] ERROR: malloc failed for slice pointers\n");
         return EXIT_FAILURE;
     }

    for (size_t k = 0; k < imageDatastr->heigth; k++)
     {
        vtkMetaInfo.dataPointer[k] = (dataType*)malloc(points_in_slice * sizeof(float));
        if (!vtkMetaInfo.dataPointer[k]) 
        {
            fprintf(stderr, "[saveNiiVolumeAsVtk] ERROR: malloc failed for slice %zu\n", k);
            for (size_t i = 0; i < k; i++) free(vtkMetaInfo.dataPointer[i]);
            free(vtkMetaInfo.dataPointer);
            return EXIT_FAILURE;
        }
     }

     //Cast NiiVolume data into dataPointer
     for (size_t k = 0; k < imageDatastr->heigth; k++) {
         void* src = imageDatastr->data[k];
         float* dst = (float*)vtkMetaInfo.dataPointer[k];

         for (size_t i = 0; i < points_in_slice; i++) {
             
             switch (imageDatastr->datatype) 
             {
             case DT_UINT8: 
             { 
                 unsigned char* p = (unsigned char*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_INT8: 
             { 
                 signed char* p = (signed char*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_UINT16: 
             { 
                 unsigned short* p = (unsigned short*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_INT16: 
             { 
                 short* p = (short*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_UINT32: 
             { 
                 unsigned int* p = (unsigned int*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_INT32: 
             { 
                 int* p = (int*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             case DT_FLOAT32: 
             { 
                 float* p = (float*)src;
                 dst[i] = p[i];
                 break;
             }
             case DT_FLOAT64: 
             { 
                 double* p = (double*)src;
                 dst[i] = (float)p[i];
                 break;
             }
             default:
                 dst[i] = 0.0f;
                 break;
             }
         }
     }

     int result = storeVtkFile(vtk_path, &vtkMetaInfo, dataForm);

     //Free the slice buffers
     for (size_t k = 0; k < imageDatastr->heigth; k++)
     {
         free(vtkMetaInfo.dataPointer[k]);
     }
     free(vtkMetaInfo.dataPointer);

     return true;
 }