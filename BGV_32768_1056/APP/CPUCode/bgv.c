#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <time.h> 
#include <sys/time.h>
#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>
#include "HErandom.h"
#include "gmp_extra_functions.h"
#include "context.h"
#include "HErandom.h"
#include "bijection_CRT_Q.h"

mpz_t p;
mpz_t *dk;
mpz_t **pk;
Context context;
int64_t timing=0;
struct timeval tvalBefore, tvalAfter;

void toRq(mpz_t data[],long length){
  for(long i=0;i<length;i++)
      mpz_mod(data[i], data[i], context.Q);
}

void toRp(mpz_t data[],long length){
  for(long i=0;i<length;i++)
      mpz_mod(data[i], data[i], p);
}

//method that returns a random sample of an element in R_q
int randRq(mpz_t result[]){
  size_t sizeBytes=mpz_sizeinbase(context.Q,8);
  for(int i=0;i<context.N;i++){
    mpz_RandomBnd(&result[i],context.Q,&sizeBytes);
    //mpz_set_str(result[i],strRandomBnd_toString(q_str),10);
  }

  return 0;
}

int randZ0n(mpz_t result[]){

  int64_t value;
  double const bignum = 0xfffffff;
  double const  bnd1=bignum/3;
  double const   bnd2= (bnd1)*2;

  for(int i=0;i<context.N;i++){
    value    = wrapRandomBnd(bignum);
    //    cout << value << endl;
    if(value<bnd1)
      mpz_set_si(result[i],-1);
    else if(value >= bnd1 && value <bnd2)
      mpz_set_si(result[i],0);
    else
      mpz_set_si(result[i],1);
  }

  return 0;
}

void sampleGaussian(mpz_t poly[])
{
  double stdev=3.2;
  double const Pi=4.0*atan(1.0); // Pi=3.1415..
  static long const bignum = 0xfffffff;

  // Use the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<context.N; i+=2) {
    double r1 = (1+wrapRandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+wrapRandomBnd(bignum))/((double)bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;
    //assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations
    // Generate two Gaussians RV's, rounded to integers
    long x = (long) floor(rr*cos(theta) +0.5);
    mpz_set_si(poly[i],x);
    if (i+1 <context.N) {
      x = (long) floor(rr*sin(theta) +0.5);
      mpz_set_si(poly[i+1],x);
    }
  }
  //poly.normalize(); //Why this?
}
/*int sampleGaussian(mpz_t result[]){
  float64_t desiredNorm = boxMuller();
  int i;
  float64_t value[N];
  for(i=0;i<context.N;i++){
    value[i]=myRandReal();
  }
  int64_t actualNorm=norm(value, N);
  int64_t returnValue[N];
 for(i=0;i<context.N;i++){
    returnValue[i]=value[i]*desiredNorm/actualNorm;
    }
 return 0;
 }*/
int keyGen(){
  mpz_t e[context.N];
  //s is dk
  //a is pk[0]
  //b is pk[1]
  initialize(e,context.N);

  randRq(pk[0]);

  sampleGaussian(dk);
  toRq(dk,context.N);//as element in Rq
  sampleGaussian(e);
  toRq(e,context.N);//as element in Rq
  CRTmultiply(pk[1],pk[0],dk,&timing);
  multiply_scalar(context.N,e,e,p);
  sum(context.N,pk[1],pk[1],e);
  toRq(pk[1],context.N);
  return 0;
}
void parseM(mpz_t m[]){
  //m \in R_p is parsed as an element in R_q with infinity norm
  //bounded by p/2
  mpz_t pHalf;
  mpz_init(pHalf);
  mpz_div_ui(pHalf,p,2);
  for (int i=0;i<context.N;i++){
    mpz_sub(m[i],m[i],pHalf);
    mpz_mmod(m[i],m[i],context.Q);
  }
}

int Enc_pk(mpz_t **c,mpz_t m[]/*,mpz_t ** r*/){
  //are c,r parameters or are them local?
  //
  //  toRq(m,N); no need to do this since p<q

  mpz_t **r;
  r=malloc(sizeof(mpz_t*)*3);
  r[0]=malloc(sizeof(mpz_t)*context.N);
  r[1]=malloc(sizeof(mpz_t)*context.N);
  r[2]=malloc(sizeof(mpz_t)*context.N);
  initialize_matrix(r,3,context.N);

  gettimeofday (&tvalBefore, NULL);
  randZ0n(r[0]);
  toRq(r[0],context.N);
  sampleGaussian(r[1]);
  toRq(r[1],context.N);
  sampleGaussian(r[2]);
  toRq(r[2],context.N);
  //in pseudocode is "b" instead of "pk[1]", but b is local to keygen
    //  parseM(m);
  gettimeofday (&tvalAfter, NULL);
  timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
   	  +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  CRTmultiply(c[0],pk[1],r[0],&timing);//,PhiPowersp,iPhiPowersp,Wpowersp,iWpowersp);
  gettimeofday (&tvalBefore, NULL);
  multiply_scalar(context.N,r[2],r[2],p);
  sum(context.N,r[2],r[2],m);
  sum(context.N,c[0],c[0],r[2]);
  toRq(c[0],context.N);
  gettimeofday (&tvalAfter, NULL);
  timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
   	  +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  CRTmultiply(c[1],r[0],pk[0],&timing);//,PhiPowersp,iPhiPowersp,Wpowersp,iWpowersp);
  gettimeofday (&tvalBefore, NULL);
  multiply_scalar(context.N,r[1],r[1],p);
  sum(context.N,c[1],c[1],r[1]);
  toRq(c[1],context.N);
  gettimeofday (&tvalAfter, NULL);
  timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
   	  +tvalAfter.tv_usec) - tvalBefore.tv_usec;

  return 0;
}

int Dec_dk(mpz_t t[],mpz_t **c){
  //s is dk
  mpz_t aux[context.N];
  initialize(aux,context.N);
  /* CRTmultiply(aux,dk,dk); */
  /* CRTmultiply(aux,aux,c[2]); */
  //c[2]==0 in these examples
  CRTmultiply(aux,dk,c[1],&timing);//,PhiPowersp,iPhiPowersp,Wpowersp,iWpowersp);
  gettimeofday (&tvalBefore, NULL);
  sub(context.N,t,c[0],aux);//c[2]==0
  toRp(t,context.N);
  gettimeofday (&tvalAfter, NULL);
  timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
     	  +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  return 0;
}

int main(){
  //printf("Importing values...\n");
  //import_values();
  import_values();
  //  test();
  //printf("Importing values OK\n");

  long long i;
  long long j;

  get_parameters(&context);

  dk=malloc(sizeof(mpz_t)*context.N);
  pk=malloc(sizeof(mpz_t*)*2);
  pk[0]=malloc(sizeof(mpz_t)*context.N);
  pk[1]=malloc(sizeof(mpz_t)*context.N);

  mpz_init(p);
  //mpz_set_str(p,  "80000000000000000003",16);// number of bits: 11
  mpz_set_str(p,"1000000000000000000003",16);
  //mpz_set_str(p,"1049089",10);// number of bits: 11
  //mpz_set_str(p,"1003",10);// number of bits: 11
  initialize(dk,context.N);
  initialize(pk[0],context.N);
  initialize(pk[1],context.N);

  //printf("Keygen...\n");
  keyGen();
  exit(0);
  //printData(dk,N);
  //printData(pk[0],N);
  //printData(pk[1],N);

  mpz_t m[context.N],m2[context.N];
  for(long i=0;i<context.N;i++){
    mpz_init_set_si(m[i],i);
    mpz_mmod(m[i],m[i],p);
    mpz_init_set_si(m2[i],context.N-i);
  }
  mpz_t **c,**c2;
  c=malloc(sizeof(mpz_t*)*3);
  c[0]=malloc(sizeof(mpz_t)*context.N);
  c[1]=malloc(sizeof(mpz_t)*context.N);
  c[2]=malloc(sizeof(mpz_t)*context.N);
  c2=malloc(sizeof(mpz_t*)*3);
  c2[0]=malloc(sizeof(mpz_t)*context.N);
  c2[1]=malloc(sizeof(mpz_t)*context.N);
  c2[2]=malloc(sizeof(mpz_t)*context.N);

  initialize_matrix(c,3,context.N);
  initialize_matrix(c2,3,context.N);
  //printf("Encrypting m...\n");
  Enc_pk(c,m);


  //printf("Encrypting m2...\n");
  Enc_pk(c2,m2);

  //printf("Adding encrypted m and m2...\n");
  gettimeofday(&tvalBefore, NULL);
  sum(context.N,c2[0],c[0],c2[0]);
  sum(context.N,c2[1],c[1],c2[1]);
  toRq(c2[0],context.N);
  toRq(c2[1],context.N);
  gettimeofday (&tvalAfter, NULL);
  timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
  	   +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  //  printData(c[0],context.N);
  mpz_t t[context.N];
  initialize(t,context.N);
  //printf("Decrypting result...\n");
  Dec_dk(t,c2);
  printf("Time in microseconds: %ld microseconds\n",
	 timing
	 );
  //printData(t,context.N);
}

