#define _USE_MATH_DEFINES

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

static int nbblocks = 1000;

static uint8_t* buffer = 0;
static double* vecl = 0;
static bool* blocks = 0;
static uint8_t* blockcol = 0;

static double* norm0 = 0;
static double* norm1 = 0;
static double* norm2 = 0;

static double* point0 = 0;
static double* point1 = 0;
static double* point2 = 0;

static int* mir = 0;

static uint8_t* stars = 0;

static int skyr = 96;
static int skyg = 149;
static int skyb = 217;

__device__ double skyfunc(double a, double b, double c, double d, double e, double f, double x)
{
	return e*x+(a*(2.0*b*x+c)*sqrt(x*(b*x+c)+d))/(4.0*b)-(a*(c*c-4.0*b*d)*log(2.0*sqrt(b)*sqrt(x*(b*x+c)+d)+2.0*b*x+c))/(8.0*sqrt(b)*sqrt(b)*sqrt(b))+f*b*x*x*x/3.0+f*c*x*x/2+f*d*x;
}

__device__ int rnbw(int nbframe)
{
	int r=0, g=0, b=0;
	double tmp;
	double x;
	
		x = fmod(nbframe * 0.006, 1.0);
		tmp = fmod(x, 1.0 / 6.0);

		if (x < 1.0/6.0)
		{
			r = 255;
			g = 1530 * tmp;
		}
		else if (x < 1.0/3.0)
		{
			g = 255;
			r = 255 - 1530 * tmp;
		}
		else if (x < 0.5)
		{
			g = 255;
			b = 1530 * tmp;
		}
		else if (x < 2.0 / 3.0)
		{
			b = 255;
			g =255 - 1530 * tmp;
		}
		else if (x < 5.0 / 6.0)
		{
			b = 255;
			r = 1530 * tmp;
		}
		else
		{
			r = 255;
			b = 255 - 1530 * tmp;
		}

		return r+256*g+256*256*b;
}

__global__ void remblock(bool* blocks,int remidx)
{
	blocks[remidx] = false;
}

__global__ void addblock(bool* blocks, int addidx)
{
	blocks[addidx] = true;
}

__global__ void changecol (uint8_t* blockcol, int buildidx, uint8_t col)
{
	blockcol[buildidx] = col;
}

__global__ void setstars(uint8_t* stars)
{
	int i;
	int tmp = blockIdx.x * blockDim.x + threadIdx.x;

	int rand = tmp;

	for (i = 0; i < 10; i++) rand = (60493 *rand+11)% 115249;

	if ((rand)%200==0)
	{
		stars[tmp] = 255*rand/ 115249;
	}
	else
	{
		stars[tmp] = 0;
	}
}

__global__ void setplanet(bool* blocks, uint8_t* blockcol)
{
	int tmp, i, j, k,l;
	int rand;


	tmp = blockIdx.x * blockDim.x + threadIdx.x;
	rand = tmp;
	i = tmp % 500;
	tmp -= i;
	tmp /= 500;
	j = tmp%1000;
	tmp -= j;
	k = tmp / 1000;

	int blockidx;
	int nbblocks = 1000;
	double disttmp;
	int tmp2;
	int tmp3;

	

	for (l = 0; l < 100; l++) rand = (60493 * rand + 11) % 479001599;

				if ((j + k) % 2 == 0) tmp2 = 2 * i;
				else tmp2 = 2 * i + 1;
				
				disttmp = sqrt((tmp2 - 500 + 0.5) * (tmp2 - 500 + 0.5) + (j - 500 + 0.5) * (j - 500 + 0.5) + (k - 500 + 0.5) * (k - 500 + 0.5));
				blockidx = i + nbblocks * j + nbblocks * nbblocks * k;



				if (disttmp < 32)
				{
					blocks[blockidx] = true;
					blockcol[blockidx] = 255;
				}
				else if (disttmp < 64)
				{
					blocks[blockidx] = true;
					blockcol[blockidx] = 215;
				}
				else if (disttmp < 120)
				{
					blocks[blockidx] = true;
					if(rand%8==0) blockcol[blockidx] = 65;
					else if(rand%8==1 || rand%8==2) blockcol[blockidx] = 29;
					else blockcol[blockidx] = 35;
				}
				else if (disttmp < 180)
				{
					blocks[blockidx] = true;
					if (rand % 6 == 0) blockcol[blockidx] = 9;
					else if (rand % 6 == 1) blockcol[blockidx] = 11;
					else if (rand % 6 == 2) blockcol[blockidx] = 16;
					else blockcol[blockidx] = 17;
				}
				else if (disttmp < 250)
				{
					blocks[blockidx] = true;
					if (rand % 6 == 0) blockcol[blockidx] = 4;
					else if (rand % 6 == 1) blockcol[blockidx] = 3;
					else blockcol[blockidx] = 5;
				}
				else if (disttmp < 497)
				{
					tmp3 = disttmp - 375;
					if (tmp3 < 0) tmp3 *= -1;
					if(tmp3>15) blocks[blockidx] = true;

					if (rand % 10 == 0) blockcol[blockidx] = 43;
					else if (rand % 10 == 1) blockcol[blockidx] = 86;
					else if (rand % 10 == 2) blockcol[blockidx] = 129;
					else if (rand % 10 == 3) blockcol[blockidx] = 172;
					else if (rand % 10 < 7) blockcol[blockidx] = 8;
					else blockcol[blockidx] = 9;

					if(rand%300==0) blockcol[blockidx] = 149;
					if (rand % 1000000 == 0)
					{
						blockcol[blockidx] = 0;
					}
				}
				else if (disttmp < 499)
				{
					blocks[blockidx] = true;
					blockcol[blockidx] = 180;

					if (k > 957) blockcol[blockidx] = 204;
					if (k > 960) blockcol[blockidx] = 215;

					if (k < 50) blockcol[blockidx] = 204;
					if (k < 48) blockcol[blockidx] = 215;
					

					if (tmp2 > 950) blockcol[blockidx] = 24;

					if (tmp2 < 75)blockcol[blockidx] = 18;

					if (j > 950) blockcol[blockidx] = 29;

					if (j < 75)blockcol[blockidx] = 101;
				}


}


__global__ void addKernel(uint8_t* buffer, double* vecl, bool* blocks, double* norm0, double* norm1, double* norm2, double* point0, double* point1, double* point2, int* mir, double pos0, double pos1, double pos2, double vec0, double vec1, double vec2, double addy0, double addy1, double addy2, double addz0, double addz1, double addz2, int nbblocks, uint8_t* blockcol, uint8_t* stars, int nbframe, int skyr,int skyg, int skyb)
{
	int i;

	double sp1 = 40000;
	double sp2 = 0.001;


	double vecn0, vecn1, vecn2;
	double cpos0, cpos1, cpos2;
	double tpos0, tpos1, tpos2;
	int cnx, cny, cnz;

	double tmin;
	double ttmp;
	double alpha;

	double px, py, pz;

	double qa, qb, qc;
	double discr;

	double t1, t2;
	double t1f, t2f;
	double tcont;

	double skyfac;

	double kb, kc, kd;

	int col, colr, colg;

	int tmp = blockIdx.x * blockDim.x + threadIdx.x;
	int tmpx = tmp % 1280;
	int tmpy = (tmp - tmpx) / 1280;

	int coll;
	int cnx2;
	int cface;

	int rnbwv;

	double u, v;

	int blockidx;
	uint8_t currblock;

	uint8_t uv;

	vecn0 = vec0 + tmpx * addy0 + tmpy * addz0;
	vecn1 = vec1 + tmpx * addy1 + tmpy * addz1;
	vecn2 = vec2 + tmpx * addy2 + tmpy * addz2;

	vecn0 /= vecl[tmp];
	vecn1 /= vecl[tmp];
	vecn2 /= vecl[tmp];


	qa = vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2;
	qb = 2 * (vecn0*pos0+ vecn1 * pos1 + vecn2 * pos2 ) - 1000 * (vecn0+vecn1+vecn2);
	qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0+pos1+pos2-500);

	discr = qb * qb - 4 * qa * qc;

	if (discr <= 0)
	{
		u = 2000*((0.5 + atan2(vecn1, vecn0) / (2.0 * M_PI)));
		v =2000*((0.5 + asin(vecn2) / M_PI));
		uv = stars[(int)u + 2000 * (int)v];

		qa = vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2;
		qb = 2 * (vecn0 * pos0 + vecn1 * pos1 + vecn2 * pos2) - 1000 * (vecn0 + vecn1 + vecn2);
		qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2) - 250000;

		discr = qb * qb - 4 * qa * qc;

		if (discr <= 0) skyfac = 0;
		else
		{
			t1 = ((-1) * qb - sqrt(discr)) / (2.0 * qa);
			t2 = ((-1) * qb + sqrt(discr)) / (2.0 * qa);

			if (t1 < 0 && t2 < 0) skyfac = 0;
			else
			{
				if (t2 > t1)
				{
					t1f = t1;
					t2f = t2;
				}
				else
				{
					t1f = t2;
					t2f = t1;
				}
				if (t1f < 0) t1f = 0;

				kb = vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2;
				kc = 2 * vecn0 * pos0 - 1000 * vecn0;
				kc += 2 * vecn1 * pos1 - 1000 * vecn1;
				kc += 2 * vecn2 * pos2 - 1000 * vecn2;
				kd = pos0 * pos0 - 1000 * (pos0-250);
				kd += pos1 * pos1 - 1000 * (pos1 - 250);
				kd += pos2 * pos2 - 1000 * (pos2 - 250);

				skyfac = skyfunc(-2000 * sp2, kb, kc, kd, 1000000*sp2, sp2, t2f) - skyfunc(-2000*sp2, kb, kc, kd, 1000000*sp2,sp2, t1f);
				skyfac /= sp1;
				if (skyfac > 1)skyfac = 1;

			}
		}

		buffer[4 * tmp] = skyfac*skyr + (1-skyfac)* uv;
		buffer[4 * tmp + 1] = skyfac * skyg + (1 - skyfac) * uv;
		buffer[4 * tmp + 2] = skyfac * skyb + (1 - skyfac) * uv;
		buffer[4 * tmp + 3] = 255;
	}
	else
	{
		t1 = ((-1) * qb - sqrt(discr)) / (2.0 * qa);
		t2 = ((-1) * qb + sqrt(discr)) / (2.0 * qa);

		if (t1 < 0 && t2 < 0)
		{
			u = 2000 * ((0.5 + atan2(vecn1, vecn0) / (2.0 * M_PI)));
			v = 2000 * ((0.5 + asin(vecn2) / M_PI));
			uv = stars[(int)u + 2000 * (int)v];

			qa = vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2;
			qb = 2 * (vecn0 * pos0 + vecn1 * pos1 + vecn2 * pos2) - 1000 * (vecn0 + vecn1 + vecn2);
			qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2) - 250000;

			discr = qb * qb - 4 * qa * qc;

			if (discr <= 0) skyfac = 0;
			else
			{
				t1 = ((-1) * qb - sqrt(discr)) / (2.0 * qa);
				t2 = ((-1) * qb + sqrt(discr)) / (2.0 * qa);

				if (t1 < 0 && t2 < 0) skyfac = 0;
				else
				{
					if (t2 > t1)
					{
						t1f = t1;
						t2f = t2;
					}
					else
					{
						t1f = t2;
						t2f = t1;
					}
					if (t1f < 0) t1f = 0;

					kb = vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2;
					kc = 2 * vecn0 * pos0 - 1000 * vecn0;
					kc += 2 * vecn1 * pos1 - 1000 * vecn1;
					kc += 2 * vecn2 * pos2 - 1000 * vecn2;
					kd = pos0 * pos0 - 1000 * (pos0 - 250);
					kd += pos1 * pos1 - 1000 * (pos1 - 250);
					kd += pos2 * pos2 - 1000 * (pos2 - 250);

					skyfac = skyfunc(-2000 * sp2, kb, kc, kd, 1000000*sp2, sp2, t2f) - skyfunc(-2000 * sp2, kb, kc, kd, 1000000*sp2, sp2, t1f);
					skyfac /= sp1;
					if (skyfac > 1)skyfac = 1;

				}
			}

			buffer[4 * tmp] = skyfac * skyr + (1 - skyfac) * uv;
			buffer[4 * tmp + 1] = skyfac * skyg + (1 - skyfac) * uv;
			buffer[4 * tmp + 2] = skyfac * skyb + (1 - skyfac) * uv;
			buffer[4 * tmp + 3] = 255;
		}
		else
		{
			if (t1 * t2 > 0)
			{
				if (t1 < t2) tcont = t1;
				else tcont = t2;

				cpos0 = pos0 + tcont * vecn0;
				cpos1 = pos1 + tcont * vecn1;
				cpos2 = pos2 + tcont * vecn2;
			}
			else
			{
				cpos0 = pos0;
				cpos1 = pos1;
				cpos2 = pos2;
			}

			px = fmod(cpos0, 1.0);
			py = fmod(cpos1, 1.0);
			pz = fmod(cpos2, 1.0);

			if (px < 0) px++;
			if (py < 0) py++;
			if (pz < 0) pz++;

			cnx = cpos0 - px;
			cny = cpos1 - py;
			cnz = cpos2 - pz;


			if ((cnx + cny + cnz) % 2 != 0)
			{
				if (px < py)
				{
					if (px < pz)
					{
						if (px < 1 - py)
						{
							if (px < 1 - pz)
							{
								px++;
								cnx--;
							}
							else
							{
								pz--;
								cnz++;
							}
						}
						else
						{
							if (1 - py < 1 - pz)
							{
								py--;
								cny++;
							}
							else
							{
								pz--;
								cnz++;
							}
						}
					}
					else
					{
						if (pz < 1 - py)
						{
							pz++;
							cnz--;
						}
						else
						{
							py--;
							cny++;
						}
					}
				}
				else
				{
					if (py < pz)
					{
						if (py < 1 - px)
						{
							if (py < 1 - pz)
							{
								py++;
								cny--;
							}
							else
							{
								pz--;
								cnz++;
							}
						}
						else
						{
							if (1 - px < 1 - pz)
							{
								px--;
								cnx++;
							}
							else
							{
								pz--;
								cnz++;
							}
						}
					}
					else
					{
						if (pz < 1 - px)
						{
							pz++;
							cnz--;
						}
						else
						{
							px--;
							cnx++;
						}
					}
				}
			}

			tpos0 = cpos0;
			tpos1 = cpos1;
			tpos2 = cpos2;

			cpos0 = px;
			cpos1 = py;
			cpos2 = pz;

			tmin = 4;


			for (i = 0; i < 12; i++)
			{
				ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
				ttmp /= norm0[i] * vecn0 + norm1[i] * vecn1 + norm2[i] * vecn2;

				if (ttmp > 0 && ttmp < tmin)
				{
					tmin = ttmp;
					coll = i;
				}
			}

			cnx -= norm0[coll];
			cny -= norm1[coll];
			cnz -= norm2[coll];

			cpos0 += tmin * vecn0;
			cpos1 += tmin * vecn1;
			cpos2 += tmin * vecn2;

			tpos0 += tmin * vecn0;
			tpos1 += tmin * vecn1;
		    tpos2 += tmin * vecn2;

			cpos0 += norm0[coll];
			cpos1 += norm1[coll];
			cpos2 += norm2[coll];

			cface = mir[coll];

			if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
			{
				currblock = 0;
			}
			else
			{
				if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
				else cnx2 = cnx / 2;

				blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

				currblock = blocks[blockidx];
			}

			while (!currblock && sqrt((tpos0-500.0)*(tpos0-500.0)+(tpos1-500.0)*(tpos1-500.0)+(tpos2-500.0)*(tpos2-500.0)) < 500.0)
			{
				tmin = 4.0;


				for (i = 0; i < 12; i++)
				{
					if (i != cface) {
						ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
						ttmp /= norm0[i] * vecn0 + norm1[i] * vecn1 + norm2[i] * vecn2;

						if (ttmp > 0 && ttmp < tmin)
						{
							tmin = ttmp;
							coll = i;
						}
					}
				}

				cnx -= norm0[coll];
				cny -= norm1[coll];
				cnz -= norm2[coll];

				cpos0 += tmin * vecn0;
				cpos1 += tmin * vecn1;
				cpos2 += tmin * vecn2;

				tpos0 += tmin * vecn0;
				tpos1 += tmin * vecn1;
				tpos2 += tmin * vecn2;

				cpos0 += norm0[coll];
				cpos1 += norm1[coll];
				cpos2 += norm2[coll];

				cface = mir[coll];

				if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
				{
					currblock = false;
				}
				else
				{
					if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
					else cnx2 = cnx / 2;
					cnx2 = cnx/2;

					blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

					currblock = blocks[blockidx];
				}
			}

			
			if (currblock)
			{
				col = blockcol[blockidx];
				if (col == 255)
				{
					rnbwv = rnbw(nbframe);

					colr = rnbwv % 256;
					rnbwv -= colr;
					rnbwv /= 256;
					colg = rnbwv % 256;
					rnbwv -= colg;
					rnbwv /= 256;
					col = rnbwv % 256;

					alpha = 1 - 0.042 * cface;
					buffer[4 * tmp] = colr * alpha;
					buffer[4 * tmp + 1] = colg * alpha;
					buffer[4 * tmp + 2] = col * alpha;
					buffer[4 * tmp + 3] = 255;
				}
				else
				{
					colr = col % 6;
					col -= colr;
					col /= 6;
					colg = col % 6;
					col -= colg;
					col /= 6;
					col %= 6;

					colr *= (255 / 5);
					colg *= (255 / 5);
					col *= (255 / 5);

					alpha = 1 - 0.042 * cface;
					buffer[4 * tmp] = colr * alpha;
					buffer[4 * tmp + 1] = colg * alpha;
					buffer[4 * tmp + 2] = col * alpha;
					buffer[4 * tmp + 3] = 255;
				}
			}
			else
			{
				buffer[4 * tmp] = skyr;
				buffer[4 * tmp + 1] = skyg;
				buffer[4 * tmp + 2] = skyb;
				buffer[4 * tmp + 3] = 255;
			}
		}
	}
}


void cudaInit(bool* blockstmp)
{
	double dist = 1;
	double sqsz = 0.01 / 4;
	int tmpx, tmpy;

	double* vecltmp = new double[1280 * 720];

	double vec0, vec1, vec2;
	double addy0, addy1, addy2;
	double addz0, addz1, addz2;
	double vecn0, vecn1, vecn2;
	double x00 = 1, x01 = 0, x02 = 0;
	double x10 = 0, x11 = 1, x12 = 0;
	double x20 = 0, x21 = 0, x22 = 1;
	double multy = (1 - 1280) * sqsz / 2;
	double multz = (720 - 1) * sqsz / 2;

	double* norm0tmp = new double[12];
	double* norm1tmp = new double[12];
	double* norm2tmp = new double[12];
	double* point0tmp = new double[12];
	double* point1tmp = new double[12];
	double* point2tmp = new double[12];

	int* mirtmp = new int[12];

	uint8_t* blockcoltmp = new uint8_t[nbblocks * nbblocks * nbblocks];
	uint8_t* starstmp = new uint8_t[2000 * 2000];


	cudaSetDevice(0);
	cudaMalloc((void**)&buffer, 4 * 1280 * 720 * sizeof(uint8_t));
	cudaMalloc((void**)&vecl, 1280 * 720 * sizeof(double));
	cudaMalloc((void**)&blocks, nbblocks * nbblocks * nbblocks * sizeof(bool));

	cudaMalloc((void**)&norm0, 12 * sizeof(double));
	cudaMalloc((void**)&norm1, 12 * sizeof(double));
	cudaMalloc((void**)&norm2, 12 * sizeof(double));

	cudaMalloc((void**)&point0, 12 * sizeof(double));
	cudaMalloc((void**)&point1, 12 * sizeof(double));
	cudaMalloc((void**)&point2, 12 * sizeof(double));

	cudaMalloc((void**)&mir, 12 * sizeof(int));

	cudaMalloc((void**)&blockcol, nbblocks*nbblocks*nbblocks * sizeof(uint8_t));

	cudaMalloc((void**)&stars, 2000*2000* sizeof(uint8_t));

	vec0 = dist * x00 + multy * x10 + multz * x20;
	vec1 = dist * x01 + multy * x11 + multz * x21;
	vec2 = dist * x02 + multy * x12 + multz * x22;

	addy0 = sqsz * x10;
	addy1 = sqsz * x11;
	addy2 = sqsz * x12;

	addz0 = -sqsz * x20;
	addz1 = -sqsz * x21;
	addz2 = -sqsz * x22;

	for (int i = 0; i < 1280 * 720; i++)
	{
		tmpx = i % 1280;
		tmpy = (i - tmpx) / 1280;

		vecn0 = vec0 + tmpx * addy0 + tmpy * addz0;
		vecn1 = vec1 + tmpx * addy1 + tmpy * addz1;
		vecn2 = vec2 + tmpx * addy2 + tmpy * addz2;

		vecltmp[i] = sqrt(vecn0 * vecn0 + vecn1 * vecn1 + vecn2 * vecn2);
	}



	point0tmp[0] = 0.5;
	point1tmp[0] = 0.5;
	point2tmp[0] = 1.5;

	point0tmp[1] = 0.5;
	point1tmp[1] = 0.5;
	point2tmp[1] = 1.5;

	point0tmp[2] = 0.5;
	point1tmp[2] = 0.5;
	point2tmp[2] = 1.5;

	point0tmp[3] = 0.5;
	point1tmp[3] = 0.5;
	point2tmp[3] = 1.5;

	point0tmp[4] = 1.5;
	point1tmp[4] = 0.5;
	point2tmp[4] = 0.5;

	point0tmp[5] = 0.5;
	point1tmp[5] = 1.5;
	point2tmp[5] = 0.5;

	point0tmp[6] = -0.5;
	point1tmp[6] = 0.5;
	point2tmp[6] = 0.5;

	point0tmp[7] = 0.5;
	point1tmp[7] = -0.5;
	point2tmp[7] = 0.5;

	point0tmp[8] = 0.5;
	point1tmp[8] = 0.5;
	point2tmp[8] = -0.5;

	point0tmp[9] = 0.5;
	point1tmp[9] = 0.5;
	point2tmp[9] = -0.5;

	point0tmp[10] = 0.5;
	point1tmp[10] = 0.5;
	point2tmp[10] = -0.5;

	point0tmp[11] = 0.5;
	point1tmp[11] = 0.5;
	point2tmp[11] = -0.5;

	norm0tmp[0] = -1;
	norm1tmp[0] = 0;
	norm2tmp[0] = -1;

	norm0tmp[1] = 0;
	norm1tmp[1] = -1;
	norm2tmp[1] = -1;

	norm0tmp[2] = 1;
	norm1tmp[2] = 0;
	norm2tmp[2] = -1;

	norm0tmp[3] = 0;
	norm1tmp[3] = 1;
	norm2tmp[3] = -1;

	norm0tmp[4] = -1;
	norm1tmp[4] = -1;
	norm2tmp[4] = 0;

	norm0tmp[5] = 1;
	norm1tmp[5] = -1;
	norm2tmp[5] = 0;

	norm0tmp[6] = 1;
	norm1tmp[6] = 1;
	norm2tmp[6] = 0;

	norm0tmp[7] = -1;
	norm1tmp[7] = 1;
	norm2tmp[7] = 0;

	norm0tmp[8] = -1;
	norm1tmp[8] = 0;
	norm2tmp[8] = 1;

	norm0tmp[9] = 0;
	norm1tmp[9] = -1;
	norm2tmp[9] = 1;

	norm0tmp[10] = 1;
	norm1tmp[10] = 0;
	norm2tmp[10] = 1;

	norm0tmp[11] = 0;
	norm1tmp[11] = 1;
	norm2tmp[11] = 1;

	mirtmp[0] = 10;
	mirtmp[1] = 11;
	mirtmp[2] = 8;
	mirtmp[3] = 9;
	mirtmp[4] = 6;
	mirtmp[5] = 7;
	mirtmp[6] = 4;
	mirtmp[7] = 5;
	mirtmp[8] = 2;
	mirtmp[9] = 3;
	mirtmp[10] = 0;
	mirtmp[11] = 1;

	cudaMemcpy(vecl, vecltmp, 1280 * 720 * sizeof(double), cudaMemcpyHostToDevice);
	
	cudaMemcpy(norm0, norm0tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(norm1, norm1tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(norm2, norm2tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(point0, point0tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(point1, point1tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(point2, point2tmp, 12 * sizeof(double), cudaMemcpyHostToDevice);


	cudaMemcpy(mir, mirtmp, 12 * sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(stars, starstmp,2000*2000 * sizeof(uint8_t), cudaMemcpyHostToDevice);

	setstars << <2000 * 2000 / 500, 500 >> > (stars);
	cudaDeviceSynchronize();

	setplanet << <nbblocks*nbblocks*nbblocks / 500, 500 >> > (blocks,blockcol);
	cudaDeviceSynchronize();

	cudaMemcpy(blockstmp, blocks, nbblocks * nbblocks * nbblocks * sizeof(bool), cudaMemcpyDeviceToHost);

}

void cudaExit()
{
	cudaFree(buffer);
	cudaFree(vecl);
	cudaDeviceReset();
}

void cudathingy(uint8_t* pixels, double pos0, double pos1, double pos2, double vec0, double vec1, double vec2, double addy0, double addy1, double addy2, double addz0, double addz1, double addz2, int remidx, int addidx, int buildidx, uint8_t col, int nbframe)
{
	if (remidx != -1)
	{
		remblock << <1, 1 >> > (blocks,remidx);
		cudaDeviceSynchronize();
	}
	if (addidx != -1)
	{
		addblock << <1, 1 >> > (blocks, addidx);
		cudaDeviceSynchronize();
		
	}
	if (buildidx != -1)
	{
		changecol << <1, 1 >> > (blockcol, buildidx,col);
		cudaDeviceSynchronize();
	}

	if (sqrt((pos0 - 500) * (pos0 - 500) + (pos1 - 500) * (pos1 - 500) + (pos2 - 500) * (pos2 - 500)) < 32)
	{
		skyr = 255;
		skyg = 174;
		skyb = 201;
	}

	addKernel <<<(int)(1280 * 720 / 600), 600 >>> (buffer, vecl, blocks, norm0, norm1, norm2, point0, point1, point2, mir, pos0, pos1, pos2, vec0, vec1, vec2, addy0, addy1, addy2, addz0, addz1, addz2, nbblocks, blockcol,stars,nbframe,skyr,skyg,skyb);
	cudaDeviceSynchronize();
	cudaMemcpy(pixels, buffer, 4 * 1280 * 720 * sizeof(uint8_t), cudaMemcpyDeviceToHost);
}