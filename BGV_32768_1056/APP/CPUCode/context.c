#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>
#include "gmp_extra_functions.h"
#include "context.h"
#include "bijection_CRT_Q.h"
#include <MaxSLiCInterface.h>
#include "Maxfiles.h"

void read_array(mpz_t *data, long long length,const  char * name);
void read_array_pointer(mpz_t **data, long long length,long long inner_length,const  char * name);
void read_array_int(int *data, long long length,const  char * name);
void read_array_int64(int64_t *data, long long length,const  char * name);

void bitReversal(mpz_t DataIn[], int64_t length);
void fillData(mpz_t dataIn1[],mpz_t dataIn2[],int64_t NHalf);
void mulByiN(int64_t length,mpz_t outputVector[length], mpz_t m,mpz_t iN);
void printPrecomputedData();
void mulByPhi(int64_t N,mpz_t dataIn1[N], mpz_t PhiPowers[N],mpz_t P);
void mulByiPhi(int64_t N,mpz_t dataIn1[N], mpz_t PhiPowers[N],mpz_t P);

int imported_values=0; //false

// Same values and implementation of the FFT for FPGA's but with everything in C. For more information about parameters check the Maxeler project.

//Wpowers are the powers of the nth primitive root of unity modulus P
int num_primes;

mpz_t *P,
  * handy_prime,
  * gen,** phi_powers,
  ** iphi_powers,
  *** W_values,
  *** iW_values,
  * iN,
  *euclidean_coeffs;

mpz_t **CRT_resultVector;
mpz_t **CRT_data1;
mpz_t **CRT_data2;
mpz_t *aux;

mpz_t Q;
int64_t N;
int n;

//only for Polynomial Multiplication using FFT-IFFT parameters
void import_precomputed_parameters(){
  //initialize data
  P=malloc(sizeof(mpz_t)*num_primes);
  iN=malloc(sizeof(mpz_t)*num_primes);
  phi_powers=(mpz_t**)malloc(sizeof(mpz_t*)*num_primes);
  iphi_powers=(mpz_t**)malloc(sizeof(mpz_t*)*num_primes);
  W_values=(mpz_t***)malloc(sizeof(mpz_t**)*num_primes);
  iW_values=(mpz_t***)malloc(sizeof(mpz_t**)*num_primes);

  int64_t i,j;
  for(i=0;i<num_primes;i++){
    phi_powers[i]=malloc(sizeof(mpz_t)*N);
    iphi_powers[i]=malloc(sizeof(mpz_t)*N);
    W_values[i]=(mpz_t**)malloc(sizeof(mpz_t*)*n);
    iW_values[i]=(mpz_t**)malloc(sizeof(mpz_t*)*n);
    initialize(phi_powers[i],N);
    initialize(iphi_powers[i],N);
    for(j=0;j<n;j++){
      W_values[i][j]=malloc(sizeof(mpz_t)*N/2);
      iW_values[i][j]=malloc(sizeof(mpz_t)*N/2);
      initialize(W_values[i][j],(N/2));
      initialize(iW_values[i][j],(N/2));
    }
  }

  initialize(P,num_primes);
  initialize(iN,num_primes);

  read_array(P,num_primes,"data/primes.dat");
  read_array(iN,num_primes,"data/iN.dat");
  char *dir= malloc(sizeof("data/data_/"));
  char *filename = malloc(sizeof("data/data_/filenameveryveryveryverylong.dat"));

  for(i=0;i<num_primes;i++){
    strcpy(dir,"data/data");
    char i_str[2];
    sprintf(i_str, "%i", (int)i);
    strcat(dir,i_str);
    strcat(dir,"/");


    //printing phi_powers
    strcpy(filename,dir);
    strcat(filename,"phi_powers.dat");
    read_array(phi_powers[i],N,filename);
    //printing iphi_powers
    strcpy(filename,dir);
    strcat(filename,"iphi_powers.dat");
    read_array(iphi_powers[i],N,filename);

    //printing W_values
    strcpy(filename,dir);
    strcat(filename,"W_values.dat");
    read_array_pointer(W_values[i],n,N/2,filename);
    //printing iW_values
    strcpy(filename,dir);
    strcat(filename,"iW_values.dat");
    read_array_pointer(iW_values[i],n,N/2,filename);
  }
}

void import_values(){
  //read number of primes

  read_array_int(&num_primes,1,"data/num_primes.dat");

  //read n and N
  int64_t nN[2];
  read_array_int64(nN,2,"data/n_N.dat");
  n=nN[0];
  N=nN[1];
  import_precomputed_parameters();
  mpz_init(Q);
  read_array(&Q,1,"data/Q.dat");
  //primes_minus_bit no necessary for C code
  //read_array(...,"data/primes_minus_bit.dat")
  //read_array(...,"data/bitsize_primes_minus_bit.dat")
  //read_array(...,"data/bitsize_primes.dat")
  euclidean_coeffs=malloc(sizeof(mpz_t)*num_primes);
  initialize(euclidean_coeffs,num_primes);
  read_array(euclidean_coeffs,num_primes,"data/coefficients_extended_euclidean.dat");


  CRT_resultVector=malloc(sizeof(mpz_t*)*num_primes);
  long i;

  for(i=0;i<num_primes;i++){
    CRT_resultVector[i]=malloc(sizeof(mpz_t)*N);
  }

  CRT_data1=malloc(sizeof(mpz_t*)*num_primes);
  for(i=0;i<num_primes;i++){
    CRT_data1[i]=malloc(sizeof(mpz_t)*N);
  }

  CRT_data2=malloc(sizeof(mpz_t*)*num_primes);
  for(i=0;i<num_primes;i++){
    CRT_data2[i]=malloc(sizeof(mpz_t)*N);
  }

  initialize_matrix(CRT_resultVector,num_primes,N);
  initialize_matrix(CRT_data1,num_primes,N);
  initialize_matrix(CRT_data2,num_primes,N);


  aux=malloc(sizeof(mpz_t)*num_primes);
  initialize(aux,num_primes);

  imported_values=1;
}

void FFT(int64_t N,int n,mpz_t  output[N],mpz_t data[N], mpz_t **W, mpz_t P){
  int64_t j,i;
  int64_t numTicks = N/2;
  mpz_t  *pTmp;
  mpz_t aux[N];
  mpz_t  *pAux=aux;
  mpz_t *pOutput=output;

  initialize(aux,N);
  if(output!=data){
    for(i=0;i<N;i++){
      mpz_set(output[i],data[i]);
    }
  }

  for(i=0;i<n;i++){
    for(j=0; j<numTicks; j++){
      mpz_mul(pAux[j],output[2*j+1],W[i][j]);
      mpz_add(pAux[j],output[2*j],pAux[j]);
      mpz_mmod(pAux[j],pAux[j],P);
      mpz_mul(pAux[j+N/2],output[2*j+1],W[i][j]);
      mpz_sub(pAux[j+N/2],output[2*j],pAux[j+N/2]);
      mpz_mmod(pAux[j+N/2],pAux[j+N/2],P);
    }
    if(i<n-1){
      pTmp = output;
      output = pAux;
      pAux = pTmp;
    }
  }

  if(n%2!=0){
    for(i = 0;i<N;i++){
      mpz_swap(pOutput[i],pAux[i]);
    }
  }else{
    for(i = 0;i<N;i++){
      mpz_swap(pOutput[i],output[i]);
    }
  }
}

void IFFT(int N,int n,mpz_t output[N],mpz_t data[N],mpz_t **iWp, mpz_t P,mpz_t iN){
  FFT(N,n,output,data,iWp,P);
  mulByiN(N,output,P,iN);
}


struct c_int128{
	uint64_t firstpart;
	int64_t secondpart;
}typedef c_int128;

void init_c_int128(c_int128* c){
  c->firstpart=0;
  c->secondpart=0;
}

struct c_uint128{
	uint64_t firstpart;
	uint64_t secondpart;
}typedef c_uint128;

void init_c_uint128(c_uint128* c){
  c->firstpart=0;
  c->secondpart=0;
}

struct c_int256{
	c_uint128 firstpart;
	c_int128 secondpart;
}typedef c_int256;

void init_c_int256(c_int256* c){
  init_c_uint128(&(c->firstpart));
  init_c_int128(&(c->secondpart));
}

struct c_uint256{
	c_uint128 firstpart;
	c_uint128 secondpart;
}typedef c_uint256;

void init_c_uint256(c_uint256* c){
  init_c_uint128(&(c->firstpart));
  init_c_uint128(&(c->secondpart));
}

struct c_int512{
	c_uint256 firstpart;
	c_int256 secondpart;
}typedef c_int512;

void init_c_int512(c_int512* c){
  init_c_uint256(&(c->firstpart));
  init_c_int256(&(c->secondpart));
}

c_int128 mpz_to_c_int128(mpz_t datum){
  c_int128 dataGo;
  init_c_int128(&dataGo);
  size_t count=2;
  mpz_export(&dataGo,&count,-1,sizeof(int64_t),0,0,datum);
  return dataGo;
}

c_uint128 mpz_to_c_uint128(mpz_t datum){
  c_uint128 dataGo;
  init_c_uint128(&dataGo);
  size_t count=2;
  mpz_export(&dataGo,&count,-1,sizeof(int64_t),0,0,datum);
  return dataGo;
}

c_int256 mpz_to_c_int256(mpz_t datum){
  c_int256 dataGo;
  init_c_int256(&dataGo);
  size_t count=4;
  mpz_export(&dataGo,&count,-1,sizeof(int64_t),0,0,datum);
  return dataGo;
}

c_uint256 mpz_to_c_uint256(mpz_t datum){
  c_uint256 dataGo;
  init_c_uint256(&dataGo);
  size_t count=4;
  mpz_export(&dataGo,&count,-1,sizeof(int64_t),0,0,datum);
  return dataGo;
}

c_int512 mpz_to_c_int512(mpz_t datum){
  c_int512 dataGo;
  init_c_int512(&dataGo);
  size_t count=8;
  mpz_export(&dataGo,&count,-1,sizeof(int64_t),0,0,datum);
  return dataGo;
}

void init_array(int64_t* dataGo,size_t bytes){
	memset(dataGo, 0, bytes);
}
void mpz_to_array(int64_t* dataGo,mpz_t datum,size_t count){//count=>number of int64_t, example if 128 bits -> count=2 (2*64=128)
	mpz_export(dataGo,&count,-1,sizeof(int64_t),0,0,datum);
}

void array_to_mpz(mpz_t output, int64_t* datum,size_t count){//count=>number of int64_t, example if 128 bits -> count=2 (2*64=128)
	mpz_import(output,count,-1,sizeof(int64_t),0,0,datum);
}

void c_int128_to_mpz(mpz_t output, c_int128 datum){
  size_t b=2;
  mpz_import(output,b,-1,sizeof(int64_t),0,0,&datum);
}

void c_uint128_to_mpz(mpz_t output,c_uint128 datum){
  size_t b=2;
  mpz_import(output,b,-1,sizeof(int64_t),0,0,&datum);
}

void c_int256_to_mpz(mpz_t output ,c_int256 datum){
  size_t b=4;
  mpz_import(output,b,-1,sizeof(int64_t),0,0,&datum);
}

void c_uint256_to_mpz(mpz_t output,c_uint256 datum){
  size_t b=4;
  mpz_import(output,b,-1,sizeof(int64_t),0,0,&datum);
}

void c_int512_to_mpz(mpz_t output ,c_int512 datum){
  size_t b=8;
  mpz_import(output,b,-1,sizeof(int64_t),0,0,&datum);
}

void printData_c_int128(c_int128 data[],int length){
	for(int i=0;i<length;i++){
		printf("data[%i]= %llx %llx\n",i,data[i].firstpart,data[i].secondpart);
	}
}

max_file_t *myMaxFile0,*myMaxFile1,*myMaxFile2,*myMaxFile3,*myMaxFile4,*myMaxFile5,*myMaxFile6,*myMaxFile7,*myMaxFile8,*myMaxFile9,*myMaxFile10,*myMaxFile11,*myMaxFile12,*myMaxFile13,*myMaxFile14,*myMaxFile15,*myMaxFile16;
int64_t context_timing=0;
size_t count=16/8; //instead of 16, it would be Poly_Mul_FFT_SizeBytes but we know that the type in DFE now is dfeInt(128)
int64_t mpz2array_timing=0,array2mpz_timing=0,bitreversal_timing;

int multiply_FPGA(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2[],int index)
{

  struct timeval tvalBefore, tvalAfter;
  int64_t context_timing2=context_timing;



  int64_t dataGo1[N * count];
  int64_t dataGo2[N * count];
  int64_t dataCome[N * count];

  init_array(dataGo1,N*count*sizeof(int64_t));
  init_array(dataGo2,N*count*sizeof(int64_t));
  init_array(dataCome,N*count*sizeof(int64_t));

  mpz_t dataFFT1[N];
  mpz_t dataFFT2[N];

  for(long i=0;i<N;i++){
    mpz_init_set(dataFFT1[i],data1[i]);
    mpz_init_set(dataFFT2[i],data2[i]);
  }
  gettimeofday (&tvalBefore, NULL);


  bitReversal(dataFFT1,N);
  bitReversal(dataFFT2,N);
  gettimeofday (&tvalAfter, NULL);

  context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
	   +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  bitreversal_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
		   +tvalAfter.tv_usec) - tvalBefore.tv_usec;

  gettimeofday (&tvalBefore, NULL);

  for (int i = 0; i < N; i++) {
    mpz_to_array(&dataGo1[i * count], dataFFT1[i], count);
    mpz_to_array(&dataGo2[i * count], dataFFT2[i], count);
  }

  gettimeofday (&tvalAfter, NULL);

  context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
	   +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  mpz2array_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
		   +tvalAfter.tv_usec) - tvalBefore.tv_usec;

//  printf("CALLING FPGA\n");
  max_engine_t *myDFE;
  long j;
  //
  // these lines are commented cause we just need to call it once in this case
  //
 if(index==0){
   if(myMaxFile0==NULL)
     myMaxFile0=Poly_Mul_FFT0_init();
   myDFE = max_load(myMaxFile0,"*");
   Poly_Mul_FFT0_actions_t actions;
   for(j=0;j<10;j++){
     printf("j:%i\n",j);
     gettimeofday (&tvalBefore, NULL);
     actions.instream_inVector1=dataGo1;
     actions.instream_inVector2=dataGo2;
     actions.outstream_resultBNH=dataCome;
     actions.instream_PhiInput=phi_input;
     actions.outstream_resultNHN=&dataCome[count*N/2];
   //  Poly_Mul_FFT0(dataGo1,dataGo2,dataCome,&dataCome[count*N/2]);
     Poly_Mul_FFT0_run(myDFE,&actions);
   }
 }
//  else if(index==1){
//    if(myMaxFile1==NULL)
//      myMaxFile1=Poly_Mul_FFT1_init();
//    myDFE = max_load(myMaxFile1,"*");
//    Poly_Mul_FFT1_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
////      actions.instream_PhiInput=phi_input;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT1_run(myDFE,&actions);
//    }
//  }else if(index==2){
//    if(myMaxFile2==NULL)
//      myMaxFile2=Poly_Mul_FFT2_init();
//    myDFE = max_load(myMaxFile2,"*");
//    Poly_Mul_FFT2_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
////      actions.instream_PhiInput=phi_input;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT2_run(myDFE,&actions);
//    }
//  }
//  else if(index==3){
//    if(myMaxFile3==NULL)
//      myMaxFile3=Poly_Mul_FFT3_init();
//    myDFE = max_load(myMaxFile3,"*");
//    Poly_Mul_FFT3_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
////      actions.instream_PhiInput=phi_input;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT3_run(myDFE,&actions);
//    }
//  }else if(index==4){
//    if(myMaxFile4==NULL)
//      myMaxFile4=Poly_Mul_FFT4_init();
//    myDFE = max_load(myMaxFile4,"*");
//    Poly_Mul_FFT4_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
////      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT4_run(myDFE,&actions);
//    }
//  }else if(index==5){
//    if(myMaxFile5==NULL)
//      myMaxFile5=Poly_Mul_FFT5_init();
//    myDFE = max_load(myMaxFile5,"*");
//    Poly_Mul_FFT5_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
////      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT5_run(myDFE,&actions);
//    }
//  }else if(index==6){
//    if(myMaxFile6==NULL)
//      myMaxFile6=Poly_Mul_FFT6_init();
//    myDFE = max_load(myMaxFile6,"*");
//    Poly_Mul_FFT6_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT6_run(myDFE,&actions);
//    }
//  }else if(index==7){
//    if(myMaxFile7==NULL)
//      myMaxFile7=Poly_Mul_FFT7_init();
//    myDFE = max_load(myMaxFile7,"*");
//    Poly_Mul_FFT7_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT7_run(myDFE,&actions);
//    }
//  }else if(index==8){
//    if(myMaxFile8==NULL)
//      myMaxFile8=Poly_Mul_FFT8_init();
//    myDFE = max_load(myMaxFile8,"*");
//    Poly_Mul_FFT8_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT8_run(myDFE,&actions);
//    }
//  }else if(index==9){
//    if(myMaxFile9==NULL)
//      myMaxFile9=Poly_Mul_FFT9_init();
//    myDFE = max_load(myMaxFile9,"*");
//    Poly_Mul_FFT9_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT9_run(myDFE,&actions);
//    }
//  }else if(index==10){
//    if(myMaxFile10==NULL)
//      myMaxFile10=Poly_Mul_FFT10_init();
//    myDFE = max_load(myMaxFile10,"*");
//    Poly_Mul_FFT10_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT10_run(myDFE,&actions);
//    }
//  }else if(index==11){
//    if(myMaxFile11==NULL)
//      myMaxFile11=Poly_Mul_FFT11_init();
//    myDFE = max_load(myMaxFile11,"*");
//    Poly_Mul_FFT11_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT11_run(myDFE,&actions);
//    }
//  }else if(index==12){
//    if(myMaxFile12==NULL)
//      myMaxFile12=Poly_Mul_FFT12_init();
//    myDFE = max_load(myMaxFile12,"*");
//    Poly_Mul_FFT12_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT12_run(myDFE,&actions);
//    }
//  }else if(index==13){
//    if(myMaxFile13==NULL)
//      myMaxFile13=Poly_Mul_FFT13_init();
//    myDFE = max_load(myMaxFile13,"*");
//    Poly_Mul_FFT13_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT13_run(myDFE,&actions);
//    }
//  }else if(index==14){
//    if(myMaxFile14==NULL)
//      myMaxFile14=Poly_Mul_FFT14_init();
//    myDFE = max_load(myMaxFile14,"*");
//    Poly_Mul_FFT14_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT14_run(myDFE,&actions);
//    }
//  }else if(index==15){
//    if(myMaxFile15==NULL)
//      myMaxFile15=Poly_Mul_FFT15_init();
//    myDFE = max_load(myMaxFile15,"*");
//    Poly_Mul_FFT15_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT15_run(myDFE,&actions);
//    }
//  }else if(index==16){
//    if(myMaxFile16==NULL)
//      myMaxFile16=Poly_Mul_FFT16_init();
//    myDFE = max_load(myMaxFile16,"*");
//    Poly_Mul_FFT16_actions_t actions;
//    for(j=0;j<10;j++){
//      gettimeofday (&tvalBefore, NULL);
//      actions.instream_inVector1=dataGo1;
//      //      actions.instream_PhiInput=phi_input;
//      actions.instream_inVector2=dataGo2;
//      actions.outstream_resultBNH=dataCome;
//      actions.outstream_resultNHN=&dataCome[count*N/2];
//      Poly_Mul_FFT16_run(myDFE,&actions);
//    }
//  }
     gettimeofday (&tvalAfter, NULL);
     max_unload(myDFE);
     printf("TIME CALLFPGA: %ld\n",((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
 			       +tvalAfter.tv_usec) - tvalBefore.tv_usec);
     context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
 	     +tvalAfter.tv_usec) - tvalBefore.tv_usec;

  //FINISHED
     //  printf("FINISHED FPGA\n");
  gettimeofday (&tvalBefore, NULL);

  for(int i=0;i<N;i++){
  	//printf("data %i,%llx %llx %llx %llx %llx %llx, ",i,dataCome[i*count],dataCome[i*count+1],dataCome[i*count+2],dataCome[i*count+3],dataCome[i*count+4],dataCome[i*count+5]);
	array_to_mpz(resultVector[i],&dataCome[i*count],count);
	//printf("dataGMP : %s\n",mpz_get_str(0,16,resultVector[i]));
  }
  gettimeofday (&tvalAfter, NULL);
//  printf("TIME ARRAY2MPZ: %ld\n",((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
//			       +tvalAfter.tv_usec) - tvalBefore.tv_usec);
  array2mpz_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
  	     +tvalAfter.tv_usec) - tvalBefore.tv_usec;
    context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
  		  +tvalAfter.tv_usec) - tvalBefore.tv_usec;

    context_timing2=context_timing-context_timing2;
  //bitReversal(resultVector,N);
  //printf("time in microseconds HERE: %lld\n",context_timing2);
  //
  //case 8192 with 225 bits minimum
  return 0;

}

int multiply(int64_t N,int64_t n,mpz_t *resultVector,mpz_t *data1,mpz_t *data2,mpz_t *phi_powers,mpz_t *iphi_powers, mpz_t **W_values,mpz_t **iW_values, mpz_t P,mpz_t iN)
{
  //some case values, one can use the precomputed_parameters.c document to generate them
  //case 8192 with 225 bits minimum
  //mpz_init_set(P,*P2);
  //mpz_init_set(iN,*iN2);
  //case 8192 with 85 bits minimum
  //mpz_init_set_str(P,"633825300114114700748356599809",10);
  //mpz_init_set_str(iN,"633747928861659364481175403935",10);
  //case 8192 with 20 bits minimum
  //mpz_init_set_str(P,"17185325057",10);
  //mpz_init_set_str(iN,"17183227239",10);
  //case 256 with 20 bits minimum
  //mpz_init_set_str(P,"536874497",10);
  //mpz_init_set_str(iN,"534777331",10);
  //case 1024 with 20 bits minimum
  //mpz_init_set_str(P,"2147756033",10);
  //mpz_init_set_str(iN,"2145658615",10);
  //case 1024 with prec_params_int64_t
  //mpz_init_set_str(P,"1054721",10);
  //mpz_init_set_str(iN,"1053691",10);
  //case 2048 with 20 bits minimum
  //mpz_init_set_str(P,"4295200769",10);
  //mpz_init_set_str(iN,"4293103503",10);
  //case 4096 with 20 bits minimum
  //mpz_init_set_str(P,"8591712257",10);
  //mpz_init_set_str(iN,"8589614671",10);


  if(!imported_values){
    //import_values();
    exit(1);
  }

  mpz_t dataFFT1[N],dataFFT2[N];
  int64_t i;
  for(i=0;i<N;i++){
    mpz_init(dataFFT1[i]);
    mpz_init(dataFFT2[i]);
    mpz_set(dataFFT1[i],data1[i]);
    mpz_set(dataFFT2[i],data2[i]);
  }

  /* printf("multiplying:\n"); */
  /*   printf("dataFFT1:\n"); */
    //printData(dataFFT1,N);
  /* printf("dataFFT2:\n"); */
  //printData(dataFFT2,N);
  //----------------------------------------------------------
  //  Polynomial multiplication using FFT over Z_p[x]/<x^n+1>
  //----------------------------------------------------------

  //first step of algorithm (multiplication of a[] and b[] by phi[])

  mulByPhi(N,dataFFT1,phi_powers,P);
  mulByPhi(N,dataFFT2,phi_powers,P);

  //FFT(dataIn1) starts
  bitReversal(dataFFT1,N);
  FFT(N,n,dataFFT1,dataFFT1,W_values,P);

  //printData(dataFFT1,N);

  //  for(i=0;i<N;i++){
  //    mpz_init_set(resultVector[i],dataFFT1[i]);
  //  }


  //FFT(dataIn1) ==> dataFFT1;

  //OBSERVATION: BOTH FFTs COULD BE EXECUTED CONCURRENTLY SINCE THEY DON'T DEPEND ON EACH OTHER

  //FFT(dataIn2) starts
  bitReversal(dataFFT2,N);

  FFT(N,n,dataFFT2,dataFFT2,W_values,P);

  //  printData(dataFFT2,N);


  //FFT(dataIn2) ==> dataFFT2;
  //FFT finished


  //FFT(dataIn1)*FFT(dataIn2)
  for(i=0;i<N;i++){
    mpz_init(resultVector[i]);
    mpz_mul(resultVector[i],(dataFFT1[i]),(dataFFT2[i]));
    mpz_mmod(resultVector[i],resultVector[i],P);
  }

  //printf("after mod:\n");
  //iFFT(actualDataOut) starts
  bitReversal(resultVector,N);
  IFFT(N,n,resultVector,resultVector,iW_values,P,iN);
  //iFFT finished

  mulByiPhi(N,resultVector,iphi_powers,P);

  //printf("result:\n");
  //printData(resultVector,N);
  return 0;
}

int multiply_CRT_values(int num_primes,int64_t N,int n,
		 mpz_t **resultVector,
		 mpz_t **data1,
		 mpz_t **data2,
		 mpz_t **phip,
		 mpz_t **iphip,
		 mpz_t ***Wp,
		 mpz_t ***iWp,
		 mpz_t *P,
		 mpz_t *iN){

  int i;
  for(i=0;i<num_primes;i++){
    //printf("iter:%i\n",i);
    //printf("data:\n");
    //printData(data1[i],N);
		  multiply(N,n,resultVector[i],data1[i],data2[i],phip[i],iphip[i],Wp[i],iWp[i],P[i],iN[i]);

    //printf("finished iter:%i\n",i);
  }

  return 0;
}

int CRTmultiply(mpz_t *resultVector,
		mpz_t *data1,
		mpz_t *data2,
		int64_t *bgv_time){

  context_timing=0;
//
  struct timeval tvalBefore, tvalAfter;


  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  long long i,j;


  gettimeofday (&tvalBefore, NULL);

  for(i=0;i<N;i++){
    from_Q_to_CRT(aux,data1[i],num_primes,P);

    for(j=0;j<num_primes;j++)
      mpz_set(CRT_data1[j][i],aux[j]);
  }

  for(i=0;i<N;i++){
    from_Q_to_CRT(aux,data2[i],num_primes,P);

    for(j=0;j<num_primes;j++)
      mpz_set(CRT_data2[j][i],aux[j]);
  }
  gettimeofday (&tvalAfter, NULL);
  printf("TIME Q2CRT:%lld\n",((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
		  +tvalAfter.tv_usec) - tvalBefore.tv_usec);
  context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
		  +tvalAfter.tv_usec) - tvalBefore.tv_usec;

  for(i=0;i<num_primes;i++){
      multiply_FPGA(N,CRT_resultVector[i],CRT_data1[i],CRT_data2[i],i);
  }

  printf("TIME MPZ2ARRAY(TOTAL): %ld\n",mpz2array_timing);
  printf("TIME ARRAY2MPZ(TOTAL): %ld\n",array2mpz_timing);
  printf("TIME BITREVERSAL(TOTAL): %ld\n",bitreversal_timing);
  array2mpz_timing=0;mpz2array_timing=0;bitreversal_timing=0;

//  for(i=0;i<num_primes;i++){
//	gettimeofday (&tvalBefore, NULL);
//    multiply(N,n,CRT_resultVector[i],CRT_data1[i],CRT_data2[i],phi_powers[i],iphi_powers[i],W_values[i],iW_values[i],P[i],iN[i]);
//    gettimeofday (&tvalAfter, NULL);
//  }



  gettimeofday (&tvalBefore, NULL);
  for(i=0;i<N;i++){

    for(j=0;j<num_primes;j++)
      mpz_set(aux[j],CRT_resultVector[j][i]);

    from_CRT_to_Q(resultVector[i],aux,num_primes,euclidean_coeffs,Q);

  }
  gettimeofday (&tvalAfter, NULL);

  context_timing+=((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
		  +tvalAfter.tv_usec) - tvalBefore.tv_usec;
  printf("TIME CRT2Q:%lld\n",((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
  		  +tvalAfter.tv_usec) - tvalBefore.tv_usec);

  printf("Time in microseconds: %lld microseconds\n",(long long int)context_timing);
(*bgv_time)+=context_timing;

  return 0;
}

int multiply_CRT(mpz_t **resultVector,
		 mpz_t **data1,
		 mpz_t **data2){

  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  return multiply_CRT_values(num_primes,N,n,resultVector,data1,data2,phi_powers,iphi_powers,W_values,iW_values,P,iN);
}

void mulByPhi(int64_t N,mpz_t dataIn1[N], mpz_t *PhiPowers,mpz_t P){
  int64_t i;
  for(i=0;i<N;i++){
    mpz_mul(dataIn1[i],dataIn1[i],PhiPowers[i]);
    mpz_mmod(dataIn1[i],dataIn1[i],P);
  }
}
void mulByiPhi(int64_t N,mpz_t dataFinal[N],mpz_t *iPhiPowers,mpz_t P){
  mulByPhi(N,dataFinal,iPhiPowers,P);
}

void mulByiN(int64_t length,mpz_t outputVector[length], mpz_t m,mpz_t iN){
  int64_t i;
  for(i=0;i<length;i++){
    mpz_mul(outputVector[i],outputVector[i],iN);
    mpz_mmod(outputVector[i],outputVector[i],m);
  }
}


void bitReversal(mpz_t x[], int64_t length){
  /* Do the bit reversal */
  int64_t i2 = length >> 1;
  mpz_t tx;
  mpz_init(tx);
  int64_t k;
  int64_t j = 0;
  int64_t i;
  for (i=0;i<length-1;i++) {
    if (i < j) {
      mpz_swap(x[i],x[j]);
    }
    k = i2;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }
}

int sum(int64_t N,mpz_t *resultVector,mpz_t *data1,mpz_t *data2){
  int64_t i;
  for(i=0;i<N;i++){
    mpz_add(resultVector[i],data1[i],data2[i]);
    //mpz_mmod(resultVector[i],resultVector[i],P);
  }
  return 0;
}

int sum_CRT(mpz_t **resultVector, mpz_t **data1, mpz_t **data2){

  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  return sum_CRT_values(num_primes,N,resultVector,data1,data2);
}
int sum_CRT_values(int num_primes, int64_t N, mpz_t **resultVector,mpz_t **data1,mpz_t **data2){
  int64_t i;
  for (i=0;i<num_primes;i++){
    sum(N,resultVector[i],data1[i],data2[i]);
  }

  return 0;
}



int sub(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2[]){
  int64_t i;
  for(i=0;i<N;i++){
    mpz_sub(resultVector[i],data1[i],data2[i]);
    //mpz_mmod(resultVector[i],resultVector[i],P);
  }
  return 0;
}

int sub_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t **data2){
  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  return sub_CRT_values(num_primes,N,resultVector,data1,data2);
}

int sub_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t **data2){
  int64_t i;
  for (i=0;i<num_primes;i++){
    sub(N,resultVector[i],data1[i],data2[i]);
  }
  return 0;
}

int sum_scalar(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2){
  int64_t i;
  for(i=0;i<N;i++){
    mpz_add(resultVector[i],data1[i],data2);
    //mpz_mmod(resultVector[i],resultVector[i],P);
  }
  return 0;
}
int sum_scalar_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t *data2){


  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  return sum_scalar_CRT_values(num_primes,N,resultVector,data1,data2);

}

int sum_scalar_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t *data2){
  int64_t i;
  for (i=0;i<num_primes;i++){
    sum_scalar(N,resultVector[i],data1[i],data2[i]);
  }
  return 0;
}

int multiply_scalar(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2){
  int64_t i;
  for(i=0;i<N;i++){
    mpz_mul(resultVector[i],data1[i],data2);
    //mpz_mmod(resultVector[i],resultVector[i],P);
  }
  return 0;
}

int multiply_scalar_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t *data2){
  int64_t i;
  for (i=0;i<num_primes;i++){
    multiply_scalar(N,resultVector[i],data1[i],data2[i]);
  }
  return 0;
}

int multiply_scalar_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t *data2){

  if(!imported_values){
    //import_values();
    exit(1);//1 means that the values have not been imported;
  }

  return multiply_scalar_CRT_values(num_primes,N,resultVector,data1,data2);

}

void read_array(mpz_t *data, long long length,const  char * name){
  FILE *file = fopen(name,"r");
  char s[1024]; //it takes numbers in hex base for upto 16384 bits
  if(file==NULL){
    fprintf(stderr,"could not open file %s\n",name);
    exit(-1);}
  long long i;
  for(i=0;i<length;i++){
    fscanf(file, "%s",s);
    mpz_set_str(data[i],s,16);
  }
  fclose(file);
}
void read_array_pointer(mpz_t **data, long long length,long long inner_length,const  char * name){
  FILE *file = fopen(name,"r");
  char s[1024]; //it takes numbers in hex base for upto 16384 bits
  if(file==NULL){
    fprintf(stderr,"could not open file %s\n",name);
    exit(-1);}
  long long i,j;
  for(i=0;i<length;i++){
    for(j=0;j<inner_length;j++){
    fscanf(file, "%s",s);
    mpz_set_str(data[i][j],s,16);
    }
  }
  fclose(file);
}
void read_array_int(int *data, long long length,const  char * name){
  FILE *file = fopen(name,"r");
  if(file==NULL){
    fprintf(stderr,"could not open file %s\n",name);
    exit(-1);}

  long long i;
  for(i=0;i<length;i++){
    fscanf(file, "%i",&data[i]);
  }
  fclose(file);
}

void read_array_int64(int64_t *data, long long length,const  char * name){
  FILE *file = fopen(name,"r");
  if(file==NULL){
    fprintf(stderr,"could not open file %s\n",name);
    exit(-1);}

  long long i;
  for(i=0;i<length;i++){
    fscanf(file, "%lli",(long long int*)&data[i]);
  }
  fclose(file);
}


void printData(mpz_t data[], int64_t length){

  char *s;
  int64_t i;
  for(i=0;i<length;i++){
    s = mpz_get_str(0,10,data[i]);
    if((i+1)%8==0) printf("\n");
    printf("data[%lld]= %s; ",(long long int)i,s);
  }
  printf("\n\n");
}

void test(){
  mpz_t **data1_CRT;
  mpz_t **data2_CRT;
  mpz_t **data3_CRT;

  data1_CRT=malloc(sizeof(mpz_t*)*num_primes);
  data2_CRT=malloc(sizeof(mpz_t*)*num_primes);
  data3_CRT=malloc(sizeof(mpz_t*)*num_primes);

  int i;
  int64_t j;
  for(i=0;i<num_primes;i++){
    data1_CRT[i]=malloc(sizeof(mpz_t)*N);
    data2_CRT[i]=malloc(sizeof(mpz_t)*N);
    data3_CRT[i]=malloc(sizeof(mpz_t)*N);
    for(j=0;j<N;j++){
      mpz_init_set_ui(data1_CRT[i][j],1);
      mpz_init_set_ui(data2_CRT[i][j],0);
      mpz_init_set_ui(data3_CRT[i][j],0);
    }
    mpz_set_ui(data2_CRT[i][N/2],1);
  }

  printf("data1_CRT:\n");
  //printData((mpz_t*)data1_CRT,num_primes*N);
  printf("data2_CRT:\n");
  //  printData((mpz_t*)data2_CRT,num_primes*N);

  multiply_CRT((mpz_t **)data3_CRT,(mpz_t**)data1_CRT,(mpz_t**)data2_CRT);


  printf("result:\n");

  mpz_t output[N];
  initialize(output,N);
  mpz_t CRT_vals[num_primes];
  mpz_init(Q);
  read_array(&Q,1,"data/Q.dat");
  initialize(CRT_vals,num_primes);

  for(j=0;j<N;j++){
    for(i=0;i<num_primes;i++){
      mpz_set(CRT_vals[i],data3_CRT[i][j]);
    }
    from_CRT_to_Q(output[j],CRT_vals,num_primes,euclidean_coeffs,Q);
  }
  printf("resultVector:\n");
  printData(output,N);
  printf("\nQ=%s\n",mpz_get_str(NULL,16,Q));
}

void get_parameters(Context* context){
  context->W_values=W_values;
  context->P=P;
  context->handy_prime=handy_prime;
  context->gen=gen;
  context->phi_powers=phi_powers;
  context->iphi_powers=iphi_powers;
  context->W_values=W_values;
  context->iW_values=iW_values;
  context->iN=iN;
  context->euclidean_coeffs=euclidean_coeffs;
  context->N=N;
  context->n=n;
  context->num_primes=num_primes;
  mpz_init_set(context->Q,Q);

}

/* 

   --------------------------------
   -----------DEPRECATED-----------
   --------------------------------

  mpz_t *PhiPowers;
  mpz_t *iPhiPowers;
  mpz_t **W;
  mpz_t **iW;

  void printData(mpz_t data[], int64_t length){

    char *s;
    for(int64_t i=0;i<length;i++){
      s = mpz_get_str(0,10,data[i]);
      if((i+1)%8==0) printf("\n");
      printf("data[%d]= %s; ",
	     i,s);
    }
    printf("\n\n");
  }
void fillData(mpz_t dataIn[],mpz_t dataIn2[], int64_t length){
  mpz_t val1,val2;
  gmp_randstate_t state,state2;
  gmp_randinit_default(state);
  gmp_randinit_default(state2);
  mpz_init(val1);
  mpz_init(val2);

  for(int64_t i=0; i < length; i++){
    mpz_urandomm(val1,state,P);
    //mpz_init_set(dataIn[i],val1);
    mpz_init_set_ui(dataIn[i],i);//rand()%Poly_Mul_FFT_P;
    mpz_urandomm(val2,state2,P);
    //mpz_init_set(dataIn2[i],val2);
    mpz_init_set_ui(dataIn2[i],0);
  }
  mpz_set_ui(dataIn2[length/2],1);

}

mpz_t resultVector[N];
mpz_t dataFFT1[N];
mpz_t dataFFT2[N];
int main()
{
  //case 8192 with 85 bits minimum
  //mpz_init_set_str(P,"633825300114114700748356599809",10);
  //mpz_init_set_str(iN,"633747928861659364481175403935",10);
  //case 8192 with 20 bits minimum
  //mpz_init_set_str(P,"17185325057",10);
  //mpz_init_set_str(iN,"17183227239",10);
  //case 256 with 20 bits minimum
  //mpz_init_set_str(P,"536874497",10);
  //mpz_init_set_str(iN,"534777331",10);
  //case 1024 with 20 bits minimum
  mpz_init_set_str(P,"2147756033",10);
  mpz_init_set_str(iN,"2145658615",10);
  //case 1024 with prec_params_int64_t
  //mpz_init_set_str(P,"1054721",10);
  //mpz_init_set_str(iN,"1053691",10);

  //case 2048 with 20 bits minimum
  //mpz_init_set_str(P,"4295200769",10);
  //mpz_init_set_str(iN,"4293103503",10);
  //case 4096 with 20 bits minimum
  //mpz_init_set_str(P,"8591712257",10);
  //mpz_init_set_str(iN,"8589614671",10);
  char *prueba;
  prueba = mpz_get_str(0,10,P);
  printf("P : %s\n",prueba);
  prueba = mpz_get_str(0,10,iN);
  printf("iN : %s\n",prueba);
  W = (mpz_t **)malloc(sizeof(mpz_t*)*n);
  iW = (mpz_t **)malloc(sizeof(mpz_t*)*n);
  getValues();

  //printPrecomputedData();

  fillData(dataFFT1,dataFFT2,N);
  struct timeval tvalBefore, tvalAfter;  // removed comma



  //----------------------------------------------------------
  //  Polynomial multiplication using FFT over Z_p[x]/<x^n+1>
  //----------------------------------------------------------


  printf("dataIn1 and dataIn2\n\n");
  printData(dataFFT1,N);
  printf("-------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
  printData(dataFFT2,N);

  gettimeofday (&tvalBefore, NULL);
  //first step of algorithm (multiplication of a[] and b[] by phi[])
  mulByPhi(dataFFT1,dataFFT2);



  //FFT(dataIn1) starts
  bitReversal(dataFFT1,N);


  Poly_Mul_FFT(N/2,dataFFT1);



  //FFT(dataIn1) ==> dataFFT1;

  //OBSERVATION: BOTH FFTs COULD BE EXECUTED CONCURRENTLY SINCE THEY DON'T DEPEND ON EACH OTHER

  //FFT(dataIn2) starts
  bitReversal(dataFFT2,N);

  Poly_Mul_FFT(N/2,dataFFT2);

  //  printData(dataFFT2,N);


  //FFT(dataIn2) ==> dataFFT2;
  //FFT finished


  //FFT(dataIn1)*FFT(dataIn2)
  for(int64_t i=0;i<N;i++){
    mpz_init(resultVector[i]);
    mpz_mul(resultVector[i],(dataFFT1[i]),(dataFFT2[i]));
    mpz_mmod(resultVector[i],resultVector[i],P);
  }
  //iFFT(actualDataOut) starts
  bitReversal(resultVector,N);

  Poly_Mul_IFFT(N/2,resultVector);

  //iFFT finished

  mulByiPhi(resultVector);

  gettimeofday (&tvalAfter, NULL);

  //algorithm finished


  printf("dataIn1*dataIn2 = \n");
  printData(resultVector,N);
  printf("Time in microseconds: %ld microseconds\n",
	 ((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000000L
	  +tvalAfter.tv_usec) - tvalBefore.tv_usec
	 );

}
*/
