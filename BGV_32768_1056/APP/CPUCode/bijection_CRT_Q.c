#include <gmp.h>

void from_CRT_to_Q(mpz_t output_number_Q,mpz_t* number_CRT,int num_primes, mpz_t *coeffs,mpz_t Q){
  int i;
  mpz_set_ui(output_number_Q,0);
  mpz_t aux;
  mpz_init(aux);
  for(i=0;i<num_primes;i++){
    mpz_mul(aux,coeffs[i],number_CRT[i]);
    mpz_add(output_number_Q,output_number_Q,aux);
  }

  mpz_mmod(output_number_Q,output_number_Q,Q);
}

void from_Q_to_CRT(mpz_t* output_number_CRT,mpz_t number_Q,int num_primes,mpz_t *P){

  int i=0;
  for(i=0;i<num_primes;i++){
    mpz_mmod(output_number_CRT[i],number_Q,P[i]);
  }
}
