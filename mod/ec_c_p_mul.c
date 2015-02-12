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
#include "bn.h"

struct point{
	u8 x[20];
	u8 y[20];
};

struct Elliptic_Curve {
	u8 p[20];
	u8 a[20];
	u8 b[20];
	struct point G;
	u8 U[20];
	u32 pad[2];
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
		if(c >= '0' && c <= '9')	t = c - '0';
		else if(c >= 'a' && c <= 'f')	t = c - 'a' + 10;
		else if(c >= 'A' && c <= 'F')	t = c - 'A' + 10;
		else				t = 0;
		res |= t << (len * 4);
	}

	return res;
}

u8 *_x_to_u8_buffer(const s8 *hex){
	u32 len = strlen(hex);
	if(len % 2 != 0) return NULL;	// (add sanity check in caller)
	
	s8 xtmp[3] = {0, 0, 0};
	u8 *res = (u8 *)malloc(sizeof(u8) * len);
	u8 *ptr = res;

	while(len--){
		xtmp[0] = *hex++;
		xtmp[1] = *hex++;
		*ptr++ = (u8) _x_to_u64(xtmp);
	}

	return res;
}

int point_is_zero(const struct point *p){
	return elt_is_zero(p->x) && elt_is_zero(p->y);
}

void point_double(struct point *r)
{	
	if(elt_is_zero(r->y)){
		bn_zero((u8 *)r, 40); return; }
		
	u8	s[20], t[20], u[20]; 

// t = px*px
	bn_mon_mul(t, r->x, r->x, EC.p, NULL);	// +aux
	
// s = 2*px*px
	bn_add(s, t, t, EC.p, 20);			// +dig, c
// s = 3*px*px
	bn_add(s, s, t, EC.p, 20);
// s = 3*px*px + a
//!!	bn_copy(tu, EC.a, 20);				//const ec_a is needed here, use (tu)
	bn_add(s, s, EC.a, EC.p, 20);		

// t = 2*py
	bn_add(t, r->y, r->y, EC.p, 20);
	
// t = 1/(2*py)
	bn_copy(u, t, 20);
	bn_mon_inv(t, u, EC.p, EC.U, NULL);	// +U[20], per_curve_constant

// s = (3*px*px+a)/(2*py)
	bn_mon_mul(s, s, t, EC.p, NULL);		// +aux
	
// rx = s*s							
	bn_copy(u, r->x, 20);				// backup old rx now ! u = rx
	bn_mon_mul(r->x, s, s, EC.p, NULL);	// +aux

// t = 2*px							reuse backed up value: u = rx
	bn_add(t, u, u, EC.p, 20);			

// rx = s*s - 2*px
	bn_sub(r->x, r->x, t, EC.p, 20);		//r->x =
	
// t = -(rx-px)						reuse backed up value: u = rx
	bn_sub(t, u, r->x, EC.p, 20);
	
// ry = -s*(rx-px)
	bn_copy(u, r->y, 20);				// backup old ry now ! u = ry
	bn_mon_mul(r->y, s, t, EC.p, NULL);	// +aux

// ry = -s*(rx-px) - py					reuse backed up value: u = ry
	bn_sub(r->y, r->y, u, EC.p, 20);		//r->y =

}	//out rx, ry

void point_add(struct point *r, const struct point *q)
{
	if(point_is_zero(r)){
		*r = *q;		//bn_copy(rx, qx, 20); bn_copy(ry, qy, 20);
		return; }

	if(point_is_zero(q)) return;
		
	u8 	s[20], t[20], u[20];

// u = qx-px
	bn_sub(u, q->x, r->x, EC.p, 20);		//u32 dig, u8 c 

	if(elt_is_zero(u)){
	// u = qy-py
		bn_sub(u, q->y, r->y, EC.p, 20);	// subs const qy !!
		
		if(elt_is_zero(u))	
			point_double(r);
		else	
			bn_zero((u8 *)r, 40);

		return;
	}

// t = 1/(qx-px)
	bn_mon_inv(t, u, EC.p, EC.U, NULL);	// +U[20], per_curve_constant
	
// u = qy-py
	bn_sub(u, q->y, r->y, EC.p, 20);		// subs const qy !!
	
// s = (qy-py)/(qx-px)
	bn_mon_mul(s, t, u, EC.p, NULL);		// +aux
	
// rx = s*s
	bn_copy(u, r->x, 20);				// backup old rx now ! u = rx
	bn_mon_mul(r->x, s, s, EC.p, NULL);	// +aux
	
// t = px+qx
	bn_add(t, u, q->x, EC.p, 20);		// adds const qx !!
// rx = s*s - (px+qx)
	bn_sub(r->x, r->x, t, EC.p, 20);

// t = -(rx-px)						reuse backed up value: u = rx
	bn_sub(t, u, r->x, EC.p, 20);
	
// ry = -s*(rx-px)
	bn_copy(u, r->y, 20);				// backup old ry now ! u = ry
	bn_mon_mul(r->y, s, t, EC.p, NULL);	// +aux
	
// ry = -s*(rx-px) - py					reuse backed up value: u = ry
	bn_sub(r->y, r->y, u, EC.p, 20);
}	//out rx, ry

/*
	let be this the main _kernel function, ideally doing P[gid] = k[gid] x EC.G 
*/
void point_mul(struct point *d, const u8 *k, const struct point *b)
{
	u8 mask;
	bn_zero((u8 *)d, 40);

	for(u8 i = 0; i < 21; i++)
		for(mask = 0x80; mask != 0; mask >>= 1){
			point_double(d);
			if((k[i] & mask) != 0)
				point_add(d, b);
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

/* precompute once u8 U[20], needed in bn_mon_inv */
void precompute(u8 *U, const u8 *N){
	u32 dig;
	u8 d[20], c;	//U has low _loc_priority
	
	bn_zero(d, 20);	d[20-1] = 2;

	c = 1; for(u8 i = 20 - 1; i < 20; i--){ dig = N[i] + 255 - d[i] + c; c = dig >> 8; U[i] = dig; }			//-	bn_sub_1(U, N, t, n);
}

int main(int argc, char *argv[])
{
	u8 *t, *Q, *k;
		
	if(!argv[1]){
	/*
		Valid test case
		----------------		
		Curve domain parameters:
		p:	c1c627e1638fdc8e24299bb041e4e23af4bb5427		is prime
		a:	c1c627e1638fdc8e24299bb041e4e23af4bb5424		divisor g = 62980
		b:	877a6d84155a1de374b72d9f9d93b36bb563b2ab		divisor g = 227169643
	*/
		t = _x_to_u8_buffer("c1c627e1638fdc8e24299bb041e4e23af4bb5427");	memcpy(EC.p, t, 20);
		t = _x_to_u8_buffer("c1c627e1638fdc8e24299bb041e4e23af4bb5424");	memcpy(EC.a, t, 20);
		t = _x_to_u8_buffer("877a6d84155a1de374b72d9f9d93b36bb563b2ab");	memcpy(EC.b, t, 20);
	/*	
		Base point:
		Gx: 010aff82b3ac72569ae645af3b527be133442131			divisor g = 32209245
		Gy: 46b8ec1e6d71e5ecb549614887d57a287df573cc			divisor g = 972
	*/
		t = _x_to_u8_buffer("010aff82b3ac72569ae645af3b527be133442131");	memcpy(EC.G.x, t, 20);
		t = _x_to_u8_buffer("46b8ec1e6d71e5ecb549614887d57a287df573cc");	memcpy(EC.G.y, t, 20);
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
		t = _x_to_u8_buffer("dfd7e09d5092e7a5d24fd2fec423f7012430ae9d");	memcpy(EC.p, t, 20);
		t = _x_to_u8_buffer("dfd7e09d5092e7a5d24fd2fec423f7012430ae9a");	memcpy(EC.a, t, 20);
		t = _x_to_u8_buffer("01914dc5f39d6da3b1fa841fdc891674fa439bd4");	memcpy(EC.b, t, 20);
	/*
		Base point:
		Gx:	70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810		divisor g = 17200
		Gy:	0b23823cd6dc3df20979373e5662f7083f6aa56f		divisor g = 30017
	*/
		t = _x_to_u8_buffer("70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810");	memcpy(EC.G.x, t, 20);
		t = _x_to_u8_buffer("0b23823cd6dc3df20979373e5662f7083f6aa56f");	memcpy(EC.G.y, t, 20);
	/*
		target T point for second curve, unknown k !!!
		T.x	5432bddd1f97418147aff016eaa6100834f2caa8		divisor g = 29512
		T.y	c498b88965689ee44df349b066cd43cbf4f2c5d0		divisor g = 3
		problem is found discrete logarithm
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
	/*
		precompute u8 U[20], per_curve_constant: bn_sub_1(U, N, t, n)
		to use in bn_mon_inv, put in _constant as EC and inv256[0x80] 	
	*/	
	precompute(t, EC.p);
	memcpy(EC.U, t, 20);
	
	// after this preparation domain_parameters remains constant!
	bn_to_mon(EC.a, EC.p, 20);
	bn_to_mon(EC.b, EC.p, 20);	//b used here
	point_to_mon(&EC.G);
/*
	printf("point:\t%db\n", sizeof (struct point));
	printf("EC:\t%db @%p %d %d\n", sizeof EC, &EC, ((long)&EC%16 == 0), ((unsigned long)&EC & 15) == 0);
	printf("inv256:\t%db @%p %d %d\n", sizeof inv256, &inv256, ((long)&inv256%16 == 0), ((unsigned long)&inv256 & 15) == 0);
*/
	struct point ec_Q;		// mon, stores pub: x+y
	memcpy(ec_Q.x, Q, 20);
	memcpy(ec_Q.y, Q + 20, 20);
	point_to_mon(&ec_Q);

	/* ec_Q prepared, can be checked against P, instead of checking after point_from_mon(&P); */
	// bn_print("ec_Q", (u8 *)&ec_Q, 40);
	
	struct point P;
	/* P = k x G */	
	point_mul(&P, k, &EC.G);
	
	bn_print("k", k, 21);		free(k);
	// bn_print("P = kG", (u8 *)&P, 40);
	
	point_from_mon(&P);
	
	bn_print("P.x", P.x, 20);
	bn_print("P.y", P.y, 20);
	
	if(memcmp((u8 *)&P, Q, 40) != 0)
		puts("FAIL!");
	else	
		puts("OK!");
	
	free(Q);
}
