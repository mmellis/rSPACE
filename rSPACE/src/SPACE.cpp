// Accumulation of all final code for wolverine simulations
// 4/22/2011 - edits since last build.

#include<fstream>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<iostream>
#include <Rcpp.h>

using namespace Rcpp;

/////////////////////////////////////////////////////////////////
//// Calculate distance between two points using latlong coordinates
/////////////////////////////////////////////////////////////////
//
//double rdist_earth(double long1, double lat1, double long2, double lat2)
//{	
//	double newdist;
//	//double PI = 3.141593;  //stored in Rcpp
//	double R = 6378.388;  //radius of earth in km
//
//	double coslat1, sinlat1, coslon1, sinlon1;
//	double coslat2, sinlat2, coslon2, sinlon2;
//
//	coslat1 = cos(lat1*PI/180); //Uses formula from R: rdist.earth()
//	sinlat1 = sin(lat1*PI/180);
//	coslon1 = cos(long1*PI/180);
//	sinlon1 = sin(long1*PI/180);
//
//	coslat2 = cos(lat2*PI/180);
//	sinlat2 = sin(lat2*PI/180);
//	coslon2 = cos(long2*PI/180);
//	sinlon2 = sin(long2*PI/180);
//
//	newdist = coslat1*coslon1*coslat2*coslon2 + 
//		coslat1*sinlon1*coslat2*sinlon2 + 
//		sinlat1*sinlat2;
//	
//	if(fabs(newdist)>1)
//		newdist = newdist/fabs(newdist);
//	newdist = R * acos(newdist);
//
//	return newdist;
//} 
//
/////////////////////////////////////////////////////////////////
//// Calculate distance between two points using km
/////////////////////////////////////////////////////////////////
//double rdist_km(double x1, double y1, double x2, double y2)
//{	
//	double newdist;
//	newdist = sqrt(pow(x2-x1,2)+pow(y2-y1,2))/1000;
//	return newdist;
//} 
//
/////////////////////////////////////////////////////////////////
//// Create a grid for the landscape
/////////////////////////////////////////////////////////////////
//
//void make_grid(double x[], double y[], double *grid_size, int *pixels, int grid[])
//{
//	double disttest = sqrt(*grid_size);
//
//	//Size of grid in km (at narrowest)
//	double dist_x, dist_y;
//	dist_x = rdist_earth(x[0], y[0], x[*pixels-1], y[0]); 
//	dist_y = rdist_earth(x[0], y[0], x[0], y[*pixels-1]);
//
//	// Size of grid in pixels (constant)
//	int n_pxls_x = 0;  //# pixels is constant, size in km changes
//	while( y[n_pxls_x] == y[0] )
//		n_pxls_x++;
//	int n_pxls_y = *pixels / n_pxls_x;
//
//	int n_cells[3];
//	n_cells[0] = n_pxls_x / int(ceil(disttest * n_pxls_x / dist_x)); //Integer division to make sure number of pixels is consistent
//	n_cells[1] = n_pxls_y / int(ceil(disttest * n_pxls_y / dist_y));
//	n_cells[2] = int(n_cells[0]*n_cells[1]); 
//
//	int origin = 0;
//	int test = 1;
//	int move = 0;
//	bool ok = 1;
//	bool test_x = 0;
//	bool test_y = 0;
//
//	for(int cn = 1; cn <= n_cells[2]; cn++)
//	{	
//		test = origin;
//		if(test < *pixels) 
//			ok = 1;
//		// Test distance from origin for each pixel
// 		while(ok)
//		{	
//			test_x = (rdist_earth(x[origin], y[test], x[test], y[test]) <= disttest);
//			test_y = (rdist_earth(x[test], y[origin], x[test], y[test]) <= disttest);
//			if(test_x && test_y)
//			{
//				grid[test] = cn;
//				test++;
//			}
//			else if(test_y)
//			{	
//				move=0;
//				while(!((x[origin] <= x[test + move]) && (y[test] > y[test + move])))
//					move++; // skips to one row down from the previous point
//				test = test + move;
//			}
//			else
//				ok = 0;
//			if(test >= *pixels)
//				ok=0;
//		}
//		// Find the next origin point 
//		if(cn % n_cells[0] > 0)
//			while(grid[origin] == cn)
//				origin++;
//		else if(cn % n_cells[0] == 0)
//		{
//			origin = test;
//			while(y[origin] == y[origin - 1])
//				origin--;
//		}
//	}
//
//	// Shift grid to center
//	for(int row = 0; row < (n_pxls_y-1); row++)
//	{
//		int ct = 0;
//		while((grid[n_pxls_x * (row + 1) - (ct + 1)] == 0) && (ct <= n_pxls_x))  //Count the # trailing zeroes
//			ct++;
//		if(ct < n_pxls_x) 
//		{ 
//			ct = ct/2;	//Integer division (5/2 = 2).
//			for(int col = n_pxls_x; col > ct; col--)
//				grid[row * n_pxls_x + col - 1] = grid[row * n_pxls_x + col - 1 - ct]; 
//			for(int col = 0; col < ct; col++)
//				grid[row * n_pxls_x + col] = 0;
//		}
//	}
//
//	*grid_size = double(n_cells[2]); // Returns the total number of cells included (0s indicate not pixels not included in any cell).
//									//Slightly larger than actual value due to rounding issues (not square).
//} 
//
/////////////////////////////////////////////////////////////////
//// Remove parts of grid that aren't in snow areas
/////////////////////////////////////////////////////////////////
//extern "C" {
//void filter_grid(int grid[], double snow[], double *cutoff, int *pixels, double *snow_cutoff)
//{
//	// Find the largest grid value (smaller than *grid_size usually)
//	int maxgrid = 0;
//	for(int i = 0; i < *pixels; i++){
//		if(grid[i] > maxgrid){
//			maxgrid = grid[i];  }}
//
//
//	//initialize n_snow & px
//	int *n_snow = new int[maxgrid + 1]; //keeps track of # pixels with snow in each grid cell
//	int *px = new int[maxgrid + 1];		//keeps track of # pixels in each grid cell
//	for(int i = 0; i <= maxgrid; i++) 
//	{
//		n_snow[i] = 0;
//		px[i] = 0;
//	}
//
//	//Count up the number of pixels in each grid cell and the number of pixels with snow.
//	for(int i = 0; i < *pixels; i++) 
//	{
//		px[grid[i]]++;
//		if(snow[i] >= *snow_cutoff) //EDITED
//			n_snow[grid[i]]++;
//	}
//	
//	// If there is too little snow in the cell, reassign grid value to 0 (= not included).
//	for(int i = 1; i <= maxgrid; i++) 
//	{
//		if(n_snow[i] <= *cutoff * px[i])
//			for(int j = 0; j < *pixels; j++)
//			{
//				if(grid[j] == i)
//					grid[j]=0;
//			}
//	}
//	delete n_snow;
//	delete px;
//}
//}
///////////////////////////////////////////////////////////////
// Sampling part - setting out homeranges on landscape
///////////////////////////////////////////////////////////////

// [[Rcpp::export]]
void sample_ind(double x[], double y[], int *N, double *buffer, int use[], int *maxid, int new_order[])
{
  int longlat;
  double (*rdist)(double, double, double, double);       
  longlat=use[new_order[0]];  // Checks whether to use rdist_earth or rdist to calculate distances
  if(longlat == 1)
    rdist = &rdist_earth;
  else
    rdist = &rdist_km;  
    
	use[new_order[0]] = 1; // Automatically include the first point in the sample = 1st wolverine
	int ct = 1; // Number of individuals included so far
	int id = 1; // Current individual to check
	double mindist;
	double newdist;

	while((ct < *N) && (id < *maxid))		// NOT <=*maxid because id starts at 0
	{ 
		mindist = 500;				// Smallest distance to any other point included in the sample
		for(int j = 0; j < id ; j++)		// Searches over all PREVIOUS points
		{
			if(use[new_order[j]])					// Calculates distance to points included in the sample already
			{
				newdist = rdist(x[new_order[j]], y[new_order[j]], x[new_order[id]], y[new_order[id]]);
				if( newdist < mindist ) // Keeps track of the smallest distance with any other point.
					mindist = newdist;
			}
		}

		if(mindist >= *buffer) //Point is only included if nearest neighbor is outside of buffer distance.
		{
			ct++;
			use[new_order[id]] = 1;
		}

		id++;
	}
	*N = ct; //Replaces # desired points with # actually fit.
}

///////////////////////////////////////////////////////////////
// Create use surface based on home range center locations
///////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void use_surface(double x_wolv[], double y_wolv[], int *N_wolv,
		 double x[], double y[], double snow[], int *pixels,
		 double *sd_x, double *sd_y, double *trunc_cutoff)
{	

	double *use = new double[*pixels]; //For the cumulative probability across all individuals.
	double *overall_use = new double[*pixels];
	const double PI = 3.141593;
	double tot;
	double term1, term2, term3;

	for(int px = 0; px < *pixels; px++)
		overall_use[px] = 1;

	for(int wolv = 0; wolv < *N_wolv; wolv++) //Builds probability of use surface for each individual.
	{
		tot = 0.0;
		for(int px = 0; px < *pixels; px++)
		{
			//Bivariate normal formula, broken up.
			term1 = pow((x[px] - x_wolv[wolv])/ *sd_x, 2) + pow((y[px] - y_wolv[wolv])/ *sd_y, 2); 
			term2 = 2* PI * *sd_x * *sd_y;
			term3 = log(term2);

			use[px] = exp((term1 + 2 * term3)/(-2));
			if(*trunc_cutoff>0)
			{
				if(term1 > pow(*trunc_cutoff,2)) //truncates home ranges
					use[px] = 0;
			}
			use[px] = snow[px] * use[px]; // rescale by snow values
			tot += use[px]; //keep track of total probability over entire surface.
		}
		for(int px = 0; px < *pixels; px++)
			overall_use[px] = overall_use[px] * (1 - use[px] / tot); //accumulates across individuals, probability of no wolverines present.
	}

	for(int px = 0; px < *pixels; px++)
		snow[px] = 1 - overall_use[px]; // 1- probability of no wolverines = probability of at least one wolverine.

	delete[] overall_use;
	delete[] use;
}

/////////////////////////////////////////////////////////////////
//// Remove individuals based on a population growth rate
/////////////////////////////////////////////////////////////////
//extern "C" {
//void reduce_pop(int IN[], int N[], double useTotal[], int *pixels, int *snowpoints,
//				double x[], double y[], double snow[], int *n_grps, double *lmda, 
//				double sd_x[], double sd_y[], double trunc_cutoff[], double *snow_cutoff)
//{
//	int snow_ct = 0;
//	int N_wolv =0;
//
//	for(int i = 1; i < *n_grps; i++)
//	{
//		N_wolv = N_wolv + N[i]; 
//	}
//	
//	int n_to_drop=0;
//	n_to_drop = (int)ceil(N_wolv * abs(1 - *lmda)-0.5);
//
//	double *x_wolv = new double[n_to_drop];
//	double *y_wolv = new double[n_to_drop];
//	double *use_tmp = new double[*pixels];
//
//	int *wolv_HRcenters = new int[N_wolv];
//	int *wolv_snowIndex = new int[N_wolv];
//	int *wolv_grp       = new int[N_wolv];
//	
//	int wolv_ct = 0;
//	for(int i = 0; i < *pixels; i++)
//	{
//		if((snow[i] >= *snow_cutoff) && (wolv_ct < N_wolv))
//		{
//			for(int j=0; j < *n_grps; j++)
//			{
//				if(IN[snow_ct + j * *snowpoints] == 1)
//				{
//					wolv_HRcenters[wolv_ct]=i; //index corresponding to x,y
//					wolv_snowIndex[wolv_ct]=snow_ct + j * *snowpoints;//index corresponding to IN
//					wolv_grp[wolv_ct] = j;
//					wolv_ct++;	
//				}
//			}
//			snow_ct++;
//		}
//		use_tmp[i] = snow[i];//Initializes a temporary holder for building use surfaces.
//	}
//
//
//	int *new_order = new int[N_wolv]; 
//	for(int i = 0; i < N_wolv; i++)
//		new_order[i] = i;
//	for(int i = 0; i < N_wolv; i++) //create random permutation of snowpoints.
//	{
//		int c = (int)((double)rand() / ((double)RAND_MAX + 1) * (N_wolv - i));
//		int t = new_order[i]; 
//		new_order[i] = new_order[i+c];
//		new_order[i+c] = t;
//	}
//	
//	for(int j=0; j < *n_grps; j++)
//	{
//		int n_drop = 0;
//		for(int k=0; k < n_to_drop; k++)
//		{
//			if(wolv_grp[new_order[k]]==j)
//			{
//				x_wolv[n_drop] = x[wolv_HRcenters[new_order[k]]];
//				y_wolv[n_drop] = y[wolv_HRcenters[new_order[k]]];
//				IN[wolv_snowIndex[new_order[k]]]=0;
//				n_drop++;
//			}
//		}
//		if(n_drop>0)
//		{
//			use_surface(x_wolv, y_wolv, &n_drop, 
//				x, y, use_tmp, pixels, &sd_x[j], &sd_y[j], &trunc_cutoff[j]); 
//			for(int i=0; i < *pixels; i++)
//			{
//				useTotal[i] = 1 - (1-useTotal[i]) / (1-use_tmp[i]);
//				use_tmp[i] = snow[i];
//			}
//		}
//		N[j] = N[j] - n_drop;
//	}
//
//	delete x_wolv;
//	delete y_wolv;
//	delete use_tmp;
//	delete wolv_HRcenters;
//	delete wolv_snowIndex;
//	delete wolv_grp;
//	delete new_order;
//}
//}
///////////////////////////////////////////////////////////////
// Decide if there were actually detections or not.
///////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void calc_prob(double use[], int grid[], int detection[], double *detectionP, int *pixels, int *max_grid, double test[])
{
	int i;
	double *det = new double[*max_grid+1]; //problem point.

	for(i = 0; i <= *max_grid; i++)
		det[i] = 1;
	for(i = 0; i < *pixels; i++)
		det[(grid[i])] = det[(grid[i])]*(1-use[i]);
	for(i = 0; i <= *max_grid; i++)
	{
		det[i] = (1-det[i])*(*detectionP);
		if(test[i] < det[i])
			detection[i] = 1;
		else
			detection[i]=0;
		test[i] = det[i]; //Returns detection probabilities for each grid index.
	}                                                                                    

	delete[] det;
}
