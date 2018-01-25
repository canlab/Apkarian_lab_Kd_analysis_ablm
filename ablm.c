#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <unistd.h>
#include <getopt.h>
#include "allocate_cao.h"
#include "nifti1.h"
#define _GNU_SOURCE

void print_usage(char *progname);
float *short2float(FILE *fp1, struct nifti_1_header *hdr_ptr);
float *new_mask(FILE *fp1, struct nifti_1_header *hdr_ptr);
float ***data_process(float *in_data, float *mask_data, struct nifti_1_header *hdr_ptr, int input_label, int radius, float thres, int thres_index, int ave_R_index, int dis_index);
float ***data_process_external_roi(float *in_data, float *mask_data, float *roi_data, struct nifti_1_header *hdr_ptr, float thres, int thres_index);
int write_image(float ***img_3D, FILE *fp2, struct nifti_1_header *hdr_ptr);

float mean(float *one_D_data, int N);
float std(float *one_D_data, int N);
float corrcoef(float *x, float *y, float x_mean, float x_std, float y_mean, float y_std, int N);
float CorrNum(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w, int radius, int input_label, int thres_index);
float CorrNum_external_roi(float *t_array, float *in_data, float *mask_data, float *roi_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w, int thres_index);
int calcu_Radius(int u, int v, int w, int radius, int i, int j, int k);
float calcu_Dis(int u, int v, int w, int i, int j, int k);
float AveR(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, int u, int v, int w);
float  Dis(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w);

int main (int argc, char **argv) 
{
  FILE *fp1, *fp2, *fp_roi, *fp3;
  struct nifti_1_header *hdr_ptr_fp1, *hdr_ptr_fp2, *hdr_ptr_roi;
  float *in_data, *mask_data, *roi_data;
  float ***img_3D;

  char *progname = argv[0];
  char *in_data_path="default", *mask_data_path="default", *roi_data_path="default", *out_data_path="default";
  float thres=0.25, radius=3;
  int input_label=0, thres_index=0;
  int c;
  int mask_index = 0;
  int external_roi_index = 0;
  int ave_R_index = 0;
  int dis_index = 0;
  
  if(argc <= 1){
	print_usage(progname);
	return(1);
  }


  while ((c = getopt (argc, argv, ":i:m:x:o:t:s:r:g:ad")) != -1)
  	switch (c)
   	{
   		case 'i':
			in_data_path  = optarg;
			break;	
		case 'm':
			mask_data_path  = optarg;
			break;	
		case 'x':
			roi_data_path  = optarg;
			external_roi_index = 1;
			break;	
		case 'o':
			out_data_path  = optarg;
			break; 
		case 't':
		        thres = atof(optarg);
			break; 
		case 's':
		        thres_index = atoi(optarg);
			break;       
		case 'r':
		        radius = atof(optarg);
			break;    
		case 'g':
		        input_label = atoi(optarg);
			break; 
		case 'a':
		        ave_R_index = 1;
			break;  
		case 'd':
		        dis_index = 1;
			break;       
   		case '?':
			if(isprint(optopt))
				fprintf (stderr, "\nError! Unknown option `-%c'.\n\n", optopt);
			else
			 	fprintf (stderr,"\nError! Unknown option character `\\x%x'.\n",optopt);
			print_usage(progname);
			exit(EXIT_FAILURE);
		case ':':
			fprintf (stderr, "\nError! Option -%c requires an argument.\n\n", optopt);
			print_usage(progname);
			exit(EXIT_FAILURE);
   		default:
			print_usage(progname);		
			exit(EXIT_FAILURE);
   	}

  if(strcmp(mask_data_path,"default")==0) {
  	mask_data_path=in_data_path;
        mask_index = 1;
  }

  if(strcmp(in_data_path,"default")==0 || strcmp(out_data_path,"default")==0){
	print_usage(progname);
	return(1);
  }
  
/* open in_data nifti file for read */
  if ( ( fp1 = fopen ( in_data_path, "rb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file %s\n", in_data_path );
    exit ( 1 );
  }

  /* read head file */
  hdr_ptr_fp1 = (struct nifti_1_header *)get_spc(1, sizeof(struct nifti_1_header));
  fread(hdr_ptr_fp1, sizeof(struct nifti_1_header), 1, fp1); 

  if ( ( fp2 = fopen ( mask_data_path, "rb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file %s\n", mask_data_path );
    exit ( 1 );
  }

  hdr_ptr_fp2 = (struct nifti_1_header *)get_spc(1, sizeof(struct nifti_1_header));
  fread(hdr_ptr_fp2, sizeof(struct nifti_1_header), 1, fp2); 
  

  /* open in_data nifti file for read again */
  rewind( fp1 );
  rewind( fp2 );

  /* open nifti file for write */
  if ( ( fp3 = fopen ( out_data_path, "wb" ) ) == NULL ) {
    fprintf ( stderr, "cannot open file %s\n", out_data_path );
    exit ( 1 );
  }

  in_data = short2float(fp1, hdr_ptr_fp1);

  if(mask_index == 0)
	mask_data = short2float(fp2, hdr_ptr_fp2);
  else
	mask_data = new_mask(fp2, hdr_ptr_fp2);

  if ( external_roi_index == 1 ) {
	if ( ( fp_roi = fopen ( roi_data_path, "rb" ) ) == NULL ) {
    		fprintf ( stderr, "cannot open file %s\n", roi_data_path );
    		exit ( 1 );
	}

	hdr_ptr_roi = (struct nifti_1_header *)get_spc(1, sizeof(struct nifti_1_header));
  	fread(hdr_ptr_roi, sizeof(struct nifti_1_header), 1, fp_roi); 
 	
	rewind( fp_roi );

  	roi_data = short2float(fp_roi, hdr_ptr_roi);
	
	img_3D=data_process_external_roi(in_data, mask_data, roi_data, hdr_ptr_fp1, thres, thres_index);

     	fclose(fp_roi);
        
  }

  else
	  img_3D=data_process(in_data, mask_data, hdr_ptr_fp1, input_label, radius, thres, thres_index, ave_R_index, dis_index);
  
  write_image(img_3D, fp3, hdr_ptr_fp1);

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);


return(0);
}
  
void print_usage(char *progname)
{
  printf("Usage: ablm  -i <input_image_name>   -m <mask_image_name>  -x <roi_image_name> -o <output_image_name>   -t <correlation threshold>  -s <correlation sign>  -r <radius>   -g <link region> -a -d\n");
  printf("     -i <input_image_name>: input_image path and name; the name is in format of .nii instead of .nii.gz\n");
  printf("     -m <mask_image_name>: calculate correlations of every voxel inside the mask with every other voxel only inside this mask; default mask image is the average input image over time; if want to use defaut mask image, ignore -m.\n");
  printf("     -x <roi_image_name>: calculate correlations of every voxel inside the roi with every other voxel only outside this roi but inside the mask.\n");
  printf("     -o <output_image_name>: output_image path and name; the name is in format of .nii instead of .nii.gz\n");
  printf("     -t <correlation threshold>: the voxel is counted only if the correlation is beyond correlation threshold; default value is 0.25 and its range is between 0 and 1.\n");
  printf("     -s <correlation sign>: 1: only count the voxel with greater than positive correlation threshold; -1: only count the voxel with less than negative correlation threshold; 0: count the voxel with above both conditions; default value is 0.\n");
  printf("     -r <sphere radius>: sphere radius and its units is voxel; default value is 3.\n");  
  printf("     -g <link region>: 1: only search the voxels outside of the sphere in mask image; -1: only search the voxels inside of the sphere in mask image ; 0: search all masked voxels; default value is 0.\n");
  printf("     -a : generate an average correlation map with correlation threshold > 0.\n");
  printf("     -d : generate an average parameter map with correlation threshold > thres .\n");
  printf("     Ex. 1. default mask image, R = 0.25, absolute R, whole brain: ablm -i filt_rsn_01.nii -o link_rsn_01.nii\n");
  printf("     Ex. 2.  Hippocampus mask image, R = 0.3, positive R, whole brain: ablm -i filt_rsn_01.nii -m Hippo_mask.nii -o link_rsn_01.nii -t 0.3 -s 1\n");
  printf("     Ex. 3.  Hippocampus mask image, R = 0.3, positive R, inside of sphere r=3 voxels link: ablm -i filt_rsn_01.nii -m Hippo_mask.nii -o link_rsn.nii -t 0.3 -s 1 -r 3 -g -1\n");
  printf("     Ex. 4.  Gray_matter mask image, Hippocampus roi image, R = 0.3, positive R: ablm -i filt_rsn_01.nii -m GM_mask.nii -x Hippo_roi.nii -o link_rsn.nii -t 0.3 -s 1\n");
  printf("     Ex. 5.  default mask image to generate an average correlation map: ablm -i filt_rsn_01.nii -o ave_map.nii -a \n");
  printf("     Ex. 6.  default mask image to generate an average distance parameter map with correlation > 0.35: ablm -i filt_rsn_01.nii -o D_map.nii -d -t 0.35\n");
   
}

float *short2float(FILE *fp1, struct nifti_1_header *hdr_ptr)
{ int i_length,j_length,k_length, t_length;
  long in_data_points;
  float *in_data;
  short *in_data_short;
  int i, bitpix;
  
  i_length = hdr_ptr->dim[1]; 
  j_length = hdr_ptr->dim[2]; 
  k_length = hdr_ptr->dim[3]; 
  t_length = hdr_ptr->dim[4]; 

  in_data_points = i_length*j_length*k_length*t_length;

  bitpix = hdr_ptr->bitpix;

  if (bitpix == 32) {
  	in_data = (float *)get_spc((size_t)(in_data_points+88), sizeof(float));
  	fread(in_data, sizeof(float), (size_t)(in_data_points+88), fp1);     
        }
  else if (bitpix == 16) {
	in_data_short = (short *)get_spc((size_t)(in_data_points+176), sizeof(short));
  	fread(in_data_short, sizeof(short), (size_t)(in_data_points+176), fp1);
	
	in_data = (float *)get_spc((size_t)(in_data_points+88), sizeof(float));
	for (i = 0; i < in_data_points; i++)
		*(in_data + i + 88) = (float)*(in_data_short+i+176);

        free(in_data_short);
     	}
  else {
    printf("The data type must be float, short or integer!\n");
    exit(1);
    }

  return (in_data);

}

float *new_mask(FILE *fp1, struct nifti_1_header *hdr_ptr)
{ int i_length,j_length,k_length, t_length;
  long temp_data_points, mask_data_points;
  float *temp_data, *mask_data, *t_array;
  int i, j, k, t;
  
  i_length = hdr_ptr->dim[1]; 
  j_length = hdr_ptr->dim[2]; 
  k_length = hdr_ptr->dim[3]; 
  t_length = hdr_ptr->dim[4]; 

  temp_data_points = i_length*j_length*k_length*t_length;
  temp_data = (float *)get_spc((size_t)(temp_data_points+88), sizeof(float));
  temp_data = short2float(fp1, hdr_ptr);

  mask_data_points = i_length*j_length*k_length;
  mask_data = (float *)get_spc((size_t)(mask_data_points+88), sizeof(float));
  t_array = (float *)get_spc((size_t)t_length, sizeof(float));
 
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	for ( t = 0; t < t_length; t++ )
     		t_array[t] = *(temp_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + 88);
   	*(mask_data+i_length*j_length*k + i_length*j + i + 88) = mean(t_array, t_length);
  }

  free(temp_data);
  free(t_array);
  return (mask_data);

}

float ***data_process(float *in_data, float *mask_data, struct nifti_1_header *hdr_ptr, int input_label, int radius, float thres, int thres_index, int ave_R_index, int dis_index)
{ int i, j, k, t;
  int i_length,j_length,k_length,t_length;
  float ***img_3D;
  float mask_voxel, *t_array, t_array_mean, t_array_std;
  int offset; 
  
  i_length = hdr_ptr->dim[1]; 
  j_length = hdr_ptr->dim[2]; 
  k_length = hdr_ptr->dim[3];
  t_length = hdr_ptr->dim[4];  


  offset = 88; /* 2816=vox_offset(bytes)*bits_per_byte=352*8, bits # for head */

  img_3D = (float ***)get_3D(i_length,j_length,k_length,sizeof(float));
  t_array = (float *)get_spc((size_t)t_length, sizeof(float));
  
/*printf("%f %d \n", thres,dis_index);*/

  /* change data to a 4-D matrix */
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
        mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);

	if(mask_voxel > 0.9) { 
		for ( t = 0; t < t_length; t++ )
     		t_array[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);
		t_array_mean = mean(t_array, t_length);
		t_array_std = std(t_array, t_length);
		
		if(ave_R_index == 0 && dis_index == 0)
			img_3D[i][j][k] = CorrNum(t_array, in_data, mask_data, t_array_mean, t_array_std, i_length, j_length, k_length, t_length, offset, thres, i, j, k, radius, input_label, thres_index);
		if(ave_R_index == 1 && dis_index == 0)
			img_3D[i][j][k] = AveR(t_array, in_data, mask_data, t_array_mean, t_array_std, i_length, j_length, k_length, t_length, offset, i, j, k);
		if(ave_R_index == 0 && dis_index == 1)
			img_3D[i][j][k] =  Dis(t_array, in_data, mask_data, t_array_mean, t_array_std, i_length, j_length, k_length, t_length, offset, thres, i, j, k);
		
        }
/*printf("%d %d %d %f\n", i,j,k,img_3D[i][j][k]);*/
  }

  free(in_data);
  free(mask_data);
  free(t_array);

  return (img_3D);
  }

float ***data_process_external_roi(float *in_data, float *mask_data, float *roi_data, struct nifti_1_header *hdr_ptr, float thres, int thres_index)
{ int i, j, k, t;
  int i_length,j_length,k_length,t_length;
  float ***img_3D;
  float mask_voxel, roi_voxel, *t_array, t_array_mean, t_array_std;
  int offset; 
  
  i_length = hdr_ptr->dim[1]; 
  j_length = hdr_ptr->dim[2]; 
  k_length = hdr_ptr->dim[3];
  t_length = hdr_ptr->dim[4];  


  offset = 88; /* 2816=vox_offset(bytes)*bits_per_byte=352*8, bits # for head */

  img_3D = (float ***)get_3D(i_length,j_length,k_length,sizeof(float));
  t_array = (float *)get_spc((size_t)t_length, sizeof(float));
  
/*printf("%f %d \n", thres,dis_index);*/

  /* change data to a 4-D matrix */
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
        mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);
	roi_voxel = *(roi_data + i_length*j_length*k + i_length*j + i + offset);

	if( roi_voxel > 0.9) { 
		for ( t = 0; t < t_length; t++ )
     			t_array[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);
		t_array_mean = mean(t_array, t_length);
		t_array_std = std(t_array, t_length);
		
		img_3D[i][j][k] = CorrNum_external_roi(t_array, in_data, mask_data, roi_data, t_array_mean, t_array_std, i_length, j_length, k_length, t_length, offset, thres, i, j, k, thres_index);

/*printf("%d %d %d %f\n", i,j,k,img_3D[i][j][k]);*/		
        }

  }

  free(in_data);
  free(mask_data);
  free(t_array);

  return (img_3D);
  }

	  
int write_image(float ***img_3D, FILE *fp2, struct nifti_1_header *hdr_ptr)
{ int i, j, k;
  int i_length,j_length,k_length;
  float *out_data;
  long out_data_points;
  float extension = 0;
  float max = 0;

  
  i_length = hdr_ptr->dim[1]; 
  j_length = hdr_ptr->dim[2]; 
  k_length = hdr_ptr->dim[3]; 

  hdr_ptr->dim[0] = 3;
  hdr_ptr->dim[4] = 1;
  hdr_ptr->datatype = 16; /*NIFTI_TYPE_FLOAT32*/
  hdr_ptr->bitpix = 32;
  
  
  out_data_points = i_length*j_length*k_length; 
  out_data = (float *)get_spc((size_t)(out_data_points), sizeof(float));
	
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
       *(out_data + i_length*j_length*k + i_length*j + i) = img_3D[i][j][k];
	if( img_3D[i][j][k] > max )
		max = img_3D[i][j][k];
  }

  hdr_ptr->cal_max = max;

  fwrite(hdr_ptr, sizeof(struct nifti_1_header), 1, fp2);

  fseek(fp2, 348, SEEK_SET);
  
  fwrite(&extension,sizeof(float),1,fp2);

  fseek(fp2, 352, SEEK_SET);

  fwrite(out_data, sizeof(float), (size_t)(out_data_points), fp2);
  
  free(out_data);
  
  return (0);
}
   
float mean(float *one_D_data, int N) {
   int i;
   float sum = 0;
   float result;
   for ( i = 0; i < N; i++ )
   	sum = sum + one_D_data[i];
   result = (float)(sum/N);
   return (result);
}

float std(float *one_D_data, int N) {
  int i;
  float one_D_mean, result, sum = 0;
  one_D_mean = mean(one_D_data, N);
  for ( i = 0; i < N; i++ )
	sum = sum + pow ((one_D_data[i] - one_D_mean), 2);
  result = sqrt (sum/(N-1));
  return(result);
}

float corrcoef(float *x, float *y, float x_mean, float x_std, float y_mean, float y_std, int N) {
  int i;
  float result, sum = 0;
  for ( i = 0; i < N; i++ )
	sum = sum + (x[i] - x_mean)*(y[i] - y_mean)/(x_std*y_std); 
  result = sum/(N-1);
  return(result);
}

float CorrNum(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w, int radius, int input_label, int thres_index) {
  int i,j,k,t;
  float *one_D_data, mask_voxel;
  float one_D_mean, one_D_std, num = 0;
  float R; 
  int ball_label;
  one_D_data = (float *)get_spc((size_t)t_length, sizeof(float));
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
	mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);
	if ( mask_voxel > 0 && ((i!=u)||(j!=v)||(k!=w)) ) {

		for ( t = 0; t < t_length; t++ )
		one_D_data[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);	
                
		one_D_mean = mean(one_D_data, t_length);
    		one_D_std = std(one_D_data, t_length);

		if(input_label == 0) {
			R = corrcoef(one_D_data, t_array, one_D_mean, one_D_std, t_array_mean, t_array_std, t_length);
			if ( thres_index == 1 && R > thres )
			num = num + 1; 
			if ( thres_index == 0 && fabs(R) > thres  )
			num = num + 1; 
			if ( thres_index == -1 && R < (-1*thres)  )
			num = num + 1; 
		}
		else {
			ball_label = calcu_Radius(u, v, w, radius, i, j, k);
			if(ball_label==input_label) {
				R = corrcoef(one_D_data, t_array, one_D_mean, one_D_std, t_array_mean, t_array_std, t_length);
				if ( thres_index == 1 && R > thres )
				num = num + 1; 
				if ( thres_index == 0 && fabs(R) > thres)
				num = num + 1; 
				if ( thres_index == -1 && R < (-1*thres) )
				num = num + 1; 
			}
		} 
	}
  }
  free(one_D_data);
  return(num);
}


float CorrNum_external_roi(float *t_array, float *in_data, float *mask_data, float *roi_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w, int thres_index) {
  int i,j,k,t;
  float *one_D_data, mask_voxel, roi_voxel;
  float one_D_mean, one_D_std, num = 0;
  float R; 
  one_D_data = (float *)get_spc((size_t)t_length, sizeof(float));
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
	mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);
	roi_voxel = *(roi_data + i_length*j_length*k + i_length*j + i + offset);
	if ( mask_voxel > 0.9 && roi_voxel == 0 && ((i!=u)||(j!=v)||(k!=w)) ) {

		for ( t = 0; t < t_length; t++ )
			one_D_data[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);	
                
		one_D_mean = mean(one_D_data, t_length);
    		one_D_std = std(one_D_data, t_length);

		R = corrcoef(one_D_data, t_array, one_D_mean, one_D_std, t_array_mean, t_array_std, t_length);
		if ( thres_index == 1 && R > thres )
			num = num + 1; 
		if ( thres_index == 0 && fabs(R) > thres  )
			num = num + 1; 
		if ( thres_index == -1 && R < (-1*thres)  )
			num = num + 1; 
	}
  }
  free(one_D_data);
  return(num);
}


float AveR(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, int u, int v, int w) {
  int i,j,k,t;
  float *one_D_data, mask_voxel, num = 0;
  float one_D_mean, one_D_std;
  float R, R_sum = 0; 

  one_D_data = (float *)get_spc((size_t)t_length, sizeof(float));
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
	mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);
	if ( mask_voxel > 0 && ((i!=u)||(j!=v)||(k!=w)) ) {

		for ( t = 0; t < t_length; t++ )
		one_D_data[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);	
                
		one_D_mean = mean(one_D_data, t_length);
    		one_D_std = std(one_D_data, t_length);
		R = corrcoef(one_D_data, t_array, one_D_mean, one_D_std, t_array_mean, t_array_std, t_length);
		if(R > 0) {
			R_sum = R_sum + R;
			num = num + 1;
		}
	}
  }
  free(one_D_data);
  return(R_sum/num);
}

float  Dis(float *t_array, float *in_data, float *mask_data, float t_array_mean, float t_array_std, int i_length, int j_length, int k_length, int t_length, int offset, float thres, int u, int v, int w) {
  int i,j,k,t;
  float *one_D_data, mask_voxel;
  float one_D_mean, one_D_std;
  float R, dis, dis_result, Dis_sum=0, num=0; 

  one_D_data = (float *)get_spc((size_t)t_length, sizeof(float));
  for ( i = 0; i < i_length; i++ )
  for ( j = 0; j < j_length; j++ ) 
  for ( k = 0; k < k_length; k++ ) {
  	
	mask_voxel = *(mask_data + i_length*j_length*k + i_length*j + i + offset);
	if ( mask_voxel > 0 && ((i!=u)||(j!=v)||(k!=w)) ) {

		for ( t = 0; t < t_length; t++ )
		one_D_data[t] = *(in_data + i_length*j_length*k_length*t + i_length*j_length*k + i_length*j + i + offset);	
                
		one_D_mean = mean(one_D_data, t_length);
    		one_D_std = std(one_D_data, t_length);
		R = corrcoef(one_D_data, t_array, one_D_mean, one_D_std, t_array_mean, t_array_std, t_length);
		if(R >= thres ) {
			dis=calcu_Dis(u, v, w, i, j, k);
			 Dis_sum =  Dis_sum + dis;
			num = num + 1;
		}
	}
  }
  if(num == 0)
	dis_result = 0;
  else
	dis_result = Dis_sum/num;

  free(one_D_data);
  return(dis_result);
}

int calcu_Radius(int u, int v, int w, int radius, int i, int j, int k) {
  int ball_label = 0;
  float temp, d1, d2, d3, distance;
   
  d1 = u - i;
  d2 = v - j;
  d3 = w - k;
  temp = pow(d1,2) + pow(d2,2) + pow(d3,2);
  distance = sqrt(temp);
  if(distance <= radius)
	ball_label = -1;	
  else  
	ball_label = 1;
  return(ball_label);
}

float calcu_Dis(int u, int v, int w, int i, int j, int k) {
  float temp, d1, d2, d3, distance;
   
  d1 = u - i;
  d2 = v - j;
  d3 = w - k;
  temp = pow(d1,2) + pow(d2,2) + pow(d3,2);
  distance = sqrt(temp);
  return(distance);
}


