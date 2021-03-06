// Copyright 2007,2008,2010  Segher Boessenkool  <segher@kernel.crashing.org>
// Licensed under the terms of the GNU GPL, version 2
// http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
/*
   further edited by masterzorag@gmail.com, 2014

   changes into buffer's managment, reducing required data to min.
   use pointers and unroll some loops to target opencl address namespace 
   ocl: we need to pass an auxiliary local_ctx across calls
   u32 dig, u8 c, d[20];
   u8 s[20], t[20], u[20];
     
   bn_add, bn_sub,		dig, c
   bn_mon_mul			dig, c, d[20]
   bn_mon_inv			dig, c, d[20], v
*/

typedef unsigned char u8;
typedef unsigned int u32;

#define LEN 20		// fix bitlen / sizeof(u8)

struct local {
	u8 d[20];
	u8 c;
	u32 dig;
};

static int bn_is_zero(const u8 *d){
	for(u8 i = 0; i < LEN; i++)
		if (d[i] != 0) return 0;

	return 1;
}

/* a _kernel user_zerofill sample */
static void bn_zero(u8 *d){
	for(u8 i = 0; i < LEN; i++) d[i] = 0;
}

/* a _kernel user_memcpy sample */
static void bn_copy(u8 *d, const u8 *a){
	for(u8 i = 0; i < LEN; i++) d[i] = a[i];
}

/* a _kernel user_memcmp sample */
static int bn_compare(const u8 *a, const u8 *b){
	for(u8 i = 0; i < LEN; i++){
		if(a[i] < b[i]) return -1;
		if(a[i] > b[i]) return 1;
	}
	return 0;
}

static void bn_print(const char *name, const u8 *a, const u32 n){
	if(name) printf("%s:\t", name);
	for(u32 i = 0; i < n; i++) printf("%02x", a[i]);
	printf("\n");
}

static void bn_add(
	u8 *d,				//io
	const u8 *a, 
	const u8 *b,			//in point_add is also qx
	const u8 *N,			//mod
	struct local *aux)		//unused
{
	u32 dig; u8 c;			//needs aux_local, aux.d[LEN] unused
	
	c = 0; for(u8 i = LEN -1; i < LEN; i--){ dig = a[i] + b[i] + c; c = dig >> 8; d[i] = dig; }				//-	c = bn_add_1(d, a, b, n);
	if(c){
		c = 1; for(u8 i = LEN -1; i < LEN; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
	}
	
	if(bn_compare(d, N) >= 0){																//-	bn_reduce(d, N, n);	
		c = 1; for(u8 i = LEN -1; i < LEN; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
	}
}

static void bn_sub(
	u8 *d,				//io
	const u8 *a,			//in point_add is also qx, qy 
	const u8 *b,
	const u8 *N,			//mod
	struct local *aux)		//unused
{
	u32 dig; u8 c;			//needs aux_local, aux.d[20] unused
	
	c = 1; for(u8 i = LEN -1; i < LEN; i--){ dig = a[i] + 255 - b[i] + c; c = dig >> 8; d[i] = dig; }		//-	c = bn_sub_1(d, a, b, n);
	c = 1 - c;
	if(c){
		c = 0; for(u8 i = LEN -1; i < LEN; i--){ dig = d[i] + N[i] + c;	c = dig >> 8; d[i] = dig; }			//-	bn_add_1(d, d, N, n);
	}
}

/* 	this lookup table must be moved in _constant address space ! */
static const u8 inv256[0x80] = {
	0x01, 0xab, 0xcd, 0xb7, 0x39, 0xa3, 0xc5, 0xef,
	0xf1, 0x1b, 0x3d, 0xa7, 0x29, 0x13, 0x35, 0xdf,
	0xe1, 0x8b, 0xad, 0x97, 0x19, 0x83, 0xa5, 0xcf,
	0xd1, 0xfb, 0x1d, 0x87, 0x09, 0xf3, 0x15, 0xbf,
	0xc1, 0x6b, 0x8d, 0x77, 0xf9, 0x63, 0x85, 0xaf,
	0xb1, 0xdb, 0xfd, 0x67, 0xe9, 0xd3, 0xf5, 0x9f,
	0xa1, 0x4b, 0x6d, 0x57, 0xd9, 0x43, 0x65, 0x8f,
	0x91, 0xbb, 0xdd, 0x47, 0xc9, 0xb3, 0xd5, 0x7f,
	0x81, 0x2b, 0x4d, 0x37, 0xb9, 0x23, 0x45, 0x6f,
	0x71, 0x9b, 0xbd, 0x27, 0xa9, 0x93, 0xb5, 0x5f,
	0x61, 0x0b, 0x2d, 0x17, 0x99, 0x03, 0x25, 0x4f,
	0x51, 0x7b, 0x9d, 0x07, 0x89, 0x73, 0x95, 0x3f,
	0x41, 0xeb, 0x0d, 0xf7, 0x79, 0xe3, 0x05, 0x2f,
	0x31, 0x5b, 0x7d, 0xe7, 0x69, 0x53, 0x75, 0x1f,
	0x21, 0xcb, 0xed, 0xd7, 0x59, 0xc3, 0xe5, 0x0f,
	0x11, 0x3b, 0x5d, 0xc7, 0x49, 0x33, 0x55, 0xff,
};

static void bn_mon_mul(u8 *io,
	const u8 *a,
	const u8 *b, 
	const u8 *N,				//mod
	struct local *aux_local		// HI _loc_priority
//	u8 d[20], c; u32 dig;		// 20 + 1 + 4 = 25b
){
	u8 d[20], c;		//needs a temp buffer !!
	u32 dig;
	
	bn_zero(d);

	for(u8 i = LEN -1; i < LEN; i--){		
		c =   -(d[LEN -1] + a[LEN -1] *b[i]) * inv256[N[LEN -1] /2];
		dig = d[LEN -1] + a[LEN -1] *b[i] + N[LEN -1] *c; dig >>= 8;
	
		for(u8 j = LEN -2; j < LEN; j--){ dig += d[j] + a[j] *b[i] + N[j] *c; d[j+1] = dig; dig >>= 8; }	
		d[0] = dig; dig >>= 8;
	
		if(dig){
			c = 1; for(u8 i = LEN -1; i < LEN; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
		}

		if(bn_compare(d, N) >= 0) {																					//-	bn_reduce(d, N, 20);
			c = 1; for(u8 i = LEN -1; i < LEN; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
		}
	}
	bn_copy(io, d);
}

static void bn_mon_inv(u8 *d,	// d = 1/a mod N
	const u8 *a,				
	const u8 *N,
	const u8 *U,				// precomputed per_curve_constant
	const u8 *V,				// precomputed per_curve_constant
	struct local *aux			// u8 d[20], c; u32 dig;
){
/*!!	now prepare d[20]: use dig, c
	bn_zero(d, 20);	d[20-1] = 1;
	
	for(u8 i = 0; i < 8 *20; i++){
		c = 0; for(u8 i = 20 - 1; i < 20; i--){ dig = d[i] + d[i] + c; c = dig >> 8; d[i] = dig; }				//-	c = bn_add_1(d, d, d, 20);
		if(c){	
			c = 1; for(u8 i = 20 - 1; i < 20; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
		}

		if(bn_compare(d, N) >= 0){
			c = 1; for(u8 i = 20 - 1; i < 20; i--){ dig = d[i] + 255 - N[i] + c; c = dig >> 8; d[i] = dig; }	//-	bn_sub_1(d, d, N, n);
		}
	}
	bn_print("d", d, 20);
*/	
	bn_copy(d, V);			//1 copy from _const to loc shall be: v = V in advance

/*	now do stuff with: d, v, use also U, a		
	as seen below, v can starts initialized per_curve_constant, 
	saving a bn_mon_mul
*/	u8 v[20];

	for(u8 i = 0; i < LEN; i++){
		for(u8 mask = 0x80; mask != 0; mask >>= 1){
			bn_mon_mul(v, d, d, N, NULL);		// +aux, v = d * d
					
			/* v can starts initialized per_curve_constant !!!
			if(mask == 0x80 && i == 0) bn_print("v", v, 20);*/
		
/* U */		if((U[i] & mask) != 0)				// per_curve_constant	
/* a */			bn_mon_mul(d, v, a, N, NULL);	// d = v * a
			else
				bn_copy(d, v);					// d = v
		}
	}
}	// out d

/* below this: not kernel functions */
void bn_to_mon(u8 *d, const u8 *N){
	for(u8 i = 0; i < 8 *LEN; i++)
		bn_add(d, d, d, N, NULL);
}

void bn_from_mon(u8 *d, const u8 *N){
	u8 t[512];

	bn_zero(t);
	t[LEN -1] = 1;
	bn_mon_mul(d, d, t, N, NULL);	//512
}