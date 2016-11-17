#ifdef __cplusplus
     extern "C" {
 #endif
       const void  mpz_RandomBnd(mpz_t* outputp, mpz_t bnd,size_t* countp);
       const char* strRandomBnd_toString(char* bnd);
       const int64_t wrapRandomBnd(double number);
 #ifdef __cplusplus
     }
 #endif
