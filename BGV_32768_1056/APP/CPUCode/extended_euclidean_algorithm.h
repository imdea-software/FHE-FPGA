// given a,b; it returns (d,x,y) such that d = x*a+y*b
// d is returned as value
int extended_euclidean_algorithm(int a,int b,int* ptr_x,int* ptr_y);
int mpz_extended_euclidean_algorithm(mpz_t mpz_a, mpz_t mpz_b, mpz_t ptr_x, mpz_t ptr_y);

