// lib/mod.c : functions

#include "hf.h" 

void ShowProgress(int i, int i_max) {
	int x, ticks=50;

	printf("\r[");
	for(x=0; x<ticks*(i/i_max); x++) printf("=");
	for(   ; x<ticks;           x++) printf(" ");
	printf("] %10d/%d", i+1, i_max);
	fflush(stdout);
}

void ReplaceStr(char *in, char *org, char *rep, char *out) {
	int n;
	char *tmp;

	tmp = strstr(in, org);
	n = tmp - in;

	strncpy(out, in, n);  
	sprintf(out+n, "%s%s", rep, tmp+strlen(org));
}

void GetDimsH5(char *fn, char *dn, hsize_t *dims) {
	hid_t file_id, dataset_id, dataspace_id;

	file_id      = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	dataset_id   = H5Dopen(file_id, dn, H5P_DEFAULT);
	dataspace_id = H5Dget_space(dataset_id);

	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

	H5Fclose(file_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
}

void ReadH5(char *fn, char *dn, double *val) {
	hid_t file_id, dataset_id;

	file_id    = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	dataset_id = H5Dopen(file_id, dn, H5P_DEFAULT);

	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(file_id);
	H5Dclose(dataset_id);
}

void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val) {
	hid_t file_id, dataset_id, dataspace_id;

	file_id      = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dataspace_id = H5Screate_simple(dim, dims, NULL);
	dataset_id   = H5Dcreate(file_id, dn, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(file_id);
	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
}
