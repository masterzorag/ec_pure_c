/*
	ec_c_p_mul.c
	2014, masterzorag@gmail.com

	entirely based on the fail0verflow C implementation of elliptic curve scalar point multiplication.
	
	this is a stripped down, pure C function to do the job, no libraries to link;
	this work leads to an opencl port, new code will follow.

	program sets curve domain parameters, an integer and a point to perform point 
	multiplication by scalar, computing the derived point.

	cryptographically speaking, it verifies private/public math correlation. 

	run without args to verify known R point for first curve:
	41da1a8f74ff8d3f1ce20ef3e9d8865c96014fe3
	73ca143c9badedf2d9d3c7573307115ccfe04f13
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bn.c"

struct point{
	u8 x[20];
	u8 y[20];
};

struct Elliptic_Curve {
	u8 p[20];
	u8 a[20];
	u8 b[20];
	struct point G;
};

static struct Elliptic_Curve EC;

typedef char s8;
typedef unsigned long long int u64;

u64 _x_to_u64(const s8 *hex){
	u64 t = 0, res = 0;
	u32 len = strlen(hex);
	char c;

	while(len--){
		c = *hex++;
		if(c >= '0' && c <= '9')
			t = c - '0';
		else if(c >= 'a' && c <= 'f')
			t = c - 'a' + 10;
		else if(c >= 'A' && c <= 'F')
			t = c - 'A' + 10;
		else
			t = 0;
		res |= t << (len * 4);
	}

	return res;
}

u8 *_x_to_u8_buffer(const s8 *hex){
	u32 len = strlen(hex);
	if(len % 2 != 0) return NULL;		//Must be aligned to 2.
	
	s8 xtmp[3] = {0, 0, 0};
	u8 *res = (u8 *)malloc(sizeof(u8) * len);
	u8 *ptr = res;

	while(len--){
		xtmp[0] = *hex++;
		xtmp[1] = *hex++;
		*ptr++ = (u8)_x_to_u64(xtmp);
	}

	return res;
}

void bn_print(const char *name, const u8 *a, const u32 n){
	printf("%s:\t", name);
	for(u32 i = 0; i < n; i++)
		printf("%02x", a[i]);

	printf("\n");
}

void elt_inv(u8 *d, u8 *a){
	u8 s[20];
	bn_copy(s, a, 20);
	bn_mon_inv(d, s, EC.p, 20);
}

int point_is_zero(struct point *p){
	return elt_is_zero(p->x) && elt_is_zero(p->y);
}

void point_double(struct point *r, const struct point *p)
{
	u8 s[20], t[20];
	struct point pp = *p;
	u8 *px, *py, *rx, *ry;

	px = pp.x;
	py = pp.y;
	rx = r->x;
	ry = r->y;

	if(elt_is_zero(py)){
		bn_zero((u8 *)r, 40);
		return;
	}

// t = px*px					elt_square(t, px);	
	bn_mon_mul(t, px, px, EC.p, 20);
// s = 2*px*px					elt_add(s, t, t);
	bn_add(s, t, t, EC.p, 20);
// s = 3*px*px					elt_add(s, s, t);
	bn_add(s, s, t, EC.p, 20);
// s = 3*px*px + a				elt_add(s, s, ec_a);	// ec_a is needed here, is const!
	bn_add(s, s, EC.a, EC.p, 20);
// t = 2*py					elt_add(t, py, py);
	bn_add(t, py, py, EC.p, 20);
// t = 1/(2*py)					
	elt_inv(t, t);
// s = (3*px*px+a)/(2*py)			elt_mul(s, s, t);
	bn_mon_mul(s, s, t, EC.p, 20);
// rx = s*s					elt_square(rx, s);
	bn_mon_mul(rx, s, s, EC.p, 20);
// t = 2*px					elt_add(t, px, px);
	bn_add(t, px, px, EC.p, 20);
// rx = s*s - 2*px				elt_sub(rx, rx, t);
	bn_sub(rx, rx, t, EC.p, 20);
	
// t = -(rx-px)					elt_sub(t, px, rx);
	bn_sub(t, px, rx, EC.p, 20);
// ry = -s*(rx-px)				elt_mul(ry, s, t);
	bn_mon_mul(ry, s, t, EC.p, 20);
// ry = -s*(rx-px) - py				elt_sub(ry, ry, py);
	bn_sub(ry, ry, py, EC.p, 20);
}

void point_add(struct point *r, const struct point *p, const struct point *q)
{
	u8 s[20], t[20], u[20];
	u8 *px, *py, *qx, *qy, *rx, *ry;
	struct point pp = *p, qq = *q;

	px = pp.x;
	py = pp.y;
	qx = qq.x;
	qy = qq.y;
	rx = r->x;
	ry = r->y;

	if(point_is_zero(&pp)){
		bn_copy(rx, qx, 20); bn_copy(ry, qy, 20); return; }

	if(point_is_zero(&qq)){
		bn_copy(rx, px, 20); bn_copy(ry, py, 20); return; }

// u = qx-px					elt_sub(u, qx, px);
	bn_sub(u, qx, px, EC.p, 20);

	if(elt_is_zero(u)){
// u = qy-py					elt_sub(u, qy, py);
		bn_sub(u, qy, py, EC.p, 20);
		
		if(elt_is_zero(u))
			point_double(r, &pp);
		else
			bn_zero((u8 *)r, 40);	//point_zero(r);

		return;
	}

// t = 1/(qx-px)				
	elt_inv(t, u);
// u = qy-py					elt_sub(u, qy, py);
	bn_sub(u, qy, py, EC.p, 20);
// s = (qy-py)/(qx-px)				elt_mul(s, t, u);
	bn_mon_mul(s, t, u, EC.p, 20);
// rx = s*s					elt_square(rx, s); -> elt_mul(d, a, a);
	bn_mon_mul(rx, s, s, EC.p, 20);
// t = px+qx					elt_add(t, px, qx);
	bn_add(t, px, qx, EC.p, 20);
// rx = s*s - (px+qx)				elt_sub(rx, rx, t);
	bn_sub(rx, rx, t, EC.p, 20);

// t = -(rx-px)					elt_sub(t, px, rx);
	bn_sub(t, px, rx, EC.p, 20);
// ry = -s*(rx-px)				elt_mul(ry, s, t);
	bn_mon_mul(ry, s, t, EC.p, 20);
// ry = -s*(rx-px) - py				elt_sub(ry, ry, py);
	bn_sub(ry, ry, py, EC.p, 20);
}	//out rx, ry

void point_mul(struct point *d, const u8 *a, const struct point *b)
{
	u32 i;
	u8 mask;

	//point_zero(d);	 //memset 0
	bn_zero((u8 *)d, 40);

	for(i = 0; i < 21; i++)
		for(mask = 0x80; mask != 0; mask >>= 1){
			point_double(d, d);
			if((a[i] & mask) != 0)
				point_add(d, d, b);
		}
}

void point_to_mon(struct point *p){
	bn_to_mon(p->x, EC.p, 20);
	bn_to_mon(p->y, EC.p, 20);
}

void point_from_mon(struct point *p){
	bn_from_mon(p->x, EC.p, 20);
	bn_from_mon(p->y, EC.p, 20);
}


int main(int argc, char *argv[])
{
	u8 *p, *a, *b, *Gx, *Gy, *Q, *k;
	
	if(!argv[1]){
	/*
		Valid test case
		----------------
		
		Curve domain parameters:
		p:	c1c627e1638fdc8e24299bb041e4e23af4bb5427		is prime
		a:	c1c627e1638fdc8e24299bb041e4e23af4bb5424		divisor g = 62980
		b:	877a6d84155a1de374b72d9f9d93b36bb563b2ab		divisor g = 227169643
	*/
		p = _x_to_u8_buffer("c1c627e1638fdc8e24299bb041e4e23af4bb5427");
		a = _x_to_u8_buffer("c1c627e1638fdc8e24299bb041e4e23af4bb5424");
		b = _x_to_u8_buffer("877a6d84155a1de374b72d9f9d93b36bb563b2ab");
	/*	
		Base point:
		Gx: 010aff82b3ac72569ae645af3b527be133442131			divisor g = 32209245
		Gy: 46b8ec1e6d71e5ecb549614887d57a287df573cc			divisor g = 972
	*/
		Gx = _x_to_u8_buffer("010aff82b3ac72569ae645af3b527be133442131");
		Gy = _x_to_u8_buffer("46b8ec1e6d71e5ecb549614887d57a287df573cc");
	/*
		known verified P point for first curve, test vectors
		P.x	41da1a8f74ff8d3f1ce20ef3e9d8865c96014fe3		divisor g = 495
		P.y	73ca143c9badedf2d9d3c7573307115ccfe04f13		divisor g = 7
	*/
		Q = _x_to_u8_buffer("41da1a8f74ff8d3f1ce20ef3e9d8865c96014fe373ca143c9badedf2d9d3c7573307115ccfe04f13");		
	/*	
		using this as 
		k: 00542d46e7b3daac8aeb81e533873aabd6d74bb710			divisor g = 37
	*/		
		k = _x_to_u8_buffer("00542d46e7b3daac8aeb81e533873aabd6d74bb710");
	} else {
	/*
		Curve domain parameters:
		p:	dfd7e09d5092e7a5d24fd2fec423f7012430ae9d		is prime
		a:	dfd7e09d5092e7a5d24fd2fec423f7012430ae9a		divisor g = 530 GCD
		b:	01914dc5f39d6da3b1fa841fdc891674fa439bd4		divisor g = 266668
		N:	00dfd7e09d5092e7a5d25167ecfcfde992ebf8ecad		is prime
	*/
		p = _x_to_u8_buffer("dfd7e09d5092e7a5d24fd2fec423f7012430ae9d");
		a = _x_to_u8_buffer("dfd7e09d5092e7a5d24fd2fec423f7012430ae9a");
		b = _x_to_u8_buffer("01914dc5f39d6da3b1fa841fdc891674fa439bd4");
	/*
		Base point:
		Gx:	70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810		divisor g = 17200
		Gy:	0b23823cd6dc3df20979373e5662f7083f6aa56f		divisor g = 30017
	*/
		Gx = _x_to_u8_buffer("70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810");
		Gy = _x_to_u8_buffer("0b23823cd6dc3df20979373e5662f7083f6aa56f");
	/*
		target T point for second curve, unknown k !!!
		T.x	5432bddd1f97418147aff016eaa6100834f2caa8		divisor g = 29512
		T.y	c498b88965689ee44df349b066cd43cbf4f2c5d0		divisor g = 3
		problem is found discrete logaritm
	*/
		Q = _x_to_u8_buffer("5432bddd1f97418147aff016eaa6100834f2caa8c498b88965689ee44df349b066cd43cbf4f2c5d0");		
	/*	
		Test: known verified P point for second curve, test vectors
		P.x	b616c81e21d66dd84906468475654cf7d6f2058a
		P.y	7338bd2600ad645b093a67f4651de9edc625295c	
		using this as 
		k: 00542d46e7b3daac8aeb81e533873aabd6d74bb710			divisor g = 37
	*/		
		k = _x_to_u8_buffer("00542d46e7b3daac8aeb81e533873aabd6d74bb710");
	}
	
	// fill global curve variables from allocated buffers and release them	
	memcpy(EC.p, p, 20);		free(p);
	memcpy(EC.a, a, 20);		free(a);
	memcpy(EC.b, b, 20);		free(b);
	memcpy(EC.G.x, Gx, 20);		free(Gx);
	memcpy(EC.G.y, Gy, 20);		free(Gy);
	
	// after this preparation domain_parameters remains constant!
	bn_to_mon(EC.a, EC.p, 20);
	bn_to_mon(EC.b, EC.p, 20);	//b used here
	point_to_mon(&EC.G);

	// known/target pub stuff: bn_print("pub", Q, 40);

	struct point ec_Q;		// mon, stores pub: x+y
	memcpy(ec_Q.x, Q, 20);
	memcpy(ec_Q.y, Q + 20, 20);
	point_to_mon(&ec_Q);

	/* ec_Q prepared, can be checked against P,
		instead of checking after point_from_mon(&P); */
	//	bn_print("ec_Q", (u8 *)&ec_Q, 40);
	
	struct point P;
	/*	P = k x G */	
	point_mul(&P, k, &EC.G);
	
	bn_print("k", k, 21);	
	free(k);
//	bn_print("P = kG", (u8 *)&P, 40);
	
	point_from_mon(&P);
	
	bn_print("P.x", (u8 *)&P.x, 20);
	bn_print("P.y", (u8 *)&P.y, 20);
	
	int res = memcmp((u8 *)&P, Q, 20);
	
	if(res != 0)
		puts("FAIL!");
	else
		printf("[%d]\n", res);
	
	free(Q);
}
