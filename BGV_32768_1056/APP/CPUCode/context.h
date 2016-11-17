void import_values();

typedef struct Context{
  int num_primes;
  mpz_t Q;
  mpz_t *P;
  mpz_t * handy_prime;
  mpz_t * gen;
  mpz_t ** phi_powers;
  mpz_t ** iphi_powers;
  mpz_t *** W_values;
  mpz_t *** iW_values;
  mpz_t * iN;
  mpz_t *euclidean_coeffs;
  int64_t N;
  int n;
}Context;

void get_parameters(Context* context);
int CRTmultiply(mpz_t *resultVector,
		mpz_t *data1,
		mpz_t *data2,
		int64_t *timing);
int multiply(int64_t N,int64_t n,mpz_t *resultVector,mpz_t *data1,mpz_t *data2,mpz_t *phi_powers,mpz_t *iphi_powers, mpz_t **W_values,mpz_t **iW_values, mpz_t P,mpz_t iN);
int multiply_scalar(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2);
int sum(int64_t N,mpz_t resultVector[N],mpz_t data1[N],mpz_t data2[N]);
int sub(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2[]);
int sum_scalar(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2);

int multiply_FPGA(int64_t N,mpz_t resultVector[],mpz_t data1[],mpz_t data2[],int index);

int multiply_CRT_values(int num_primes,int64_t N,int n,
		 mpz_t **resultVector,
		 mpz_t **data1,
		 mpz_t **data2,
		 mpz_t **phip,
		 mpz_t **iphip,
		 mpz_t ***Wp,
		 mpz_t ***iWp,
		 mpz_t *P,
		 mpz_t *iN);

int multiply_CRT(mpz_t **resultVector,
		 mpz_t **data1,
		 mpz_t **data2);

int sum_scalar_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t *data2);
int sum_scalar_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t *data2);
int multiply_scalar_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t *data2);
int multiply_scalar_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t *data2);
int sum_CRT_values(int num_primes, int64_t N, mpz_t **resultVector,mpz_t **data1,mpz_t **data2);
int sum_CRT(mpz_t **resultVector, mpz_t **data1, mpz_t **data2);
int sub_CRT_values(int num_primes, int64_t N,mpz_t **resultVector,mpz_t **data1,mpz_t **data2);
int sub_CRT(mpz_t **resultVector,mpz_t **data1,mpz_t **data2);
void test();
