#include <gmp.h>
#include <inttypes.h>

void initialize(mpz_t* vector,int64_t length){
  int64_t i;
  for(i=0;i<length;i++){
    mpz_init(vector[i]);
  }
}

void initialize_matrix(mpz_t **vector, int64_t outerlength, int64_t innerlength){
  int64_t i;
  for(i=0;i<outerlength;i++)
    initialize(vector[i],innerlength);
}
