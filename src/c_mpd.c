#include <stdlib.h>
#include <R.h>

unsigned int get_2d (unsigned int a, unsigned int b, unsigned int k)
{
	return a + b * k;
}

/*
	compute community distance
	this is the entry point from R
	
	species, sites: paired array, listing site species pairs
		should be indexed to 0
		so sites[0] = 1, species[0] = 156 means sp 156 is present at site 0
	N: length of species and site arrays
	dist: 2-d array of distances between species (column major order; use get_2d to index)
	K: number of species
	dest: destination array, of length S*S
	S: number of sites
*/
void c_mpd (int * species, int * sites, int * N, double * dist, int * K, double * dest, int * S)
{
	unsigned int n = *N;
	unsigned int k = *K;
	unsigned int s = *S;
	
	// create an array to keep track of the number of species compared at each site pair
	unsigned int * ncomp = malloc( (s*s) * sizeof *ncomp);
	for(int i = 0; i < s*s; ++i) ncomp[i] = 0;
	
	for(int i = 0; i < n; ++i)
	{
		unsigned int si1 = sites[i];
		unsigned int sp1 = species[i];
		for(int j = i; j < n; ++j)
		{
			unsigned int si2 = sites[j];
			unsigned int sp2 = species[j];
			double d = dist[get_2d(sp1,sp2,k)];
			dest[get_2d(si1, si2, s)] += d;
			(ncomp[get_2d(si1, si2, s)])++;
			if(si1 != si2)
			{
				dest[get_2d(si2, si1, s)] += d;
				(ncomp[get_2d(si2, si1, s)])++;
			}
		}
	}
	
	for(int i = 0; i < s*s; ++i)
		dest[i] /= ncomp[i];
	free(ncomp);
}
