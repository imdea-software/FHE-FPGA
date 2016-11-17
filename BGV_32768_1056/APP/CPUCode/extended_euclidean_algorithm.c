#include <gmp.h>

/* int main(void){ */

/*   int a=7,b=5; */
/*   int x,y; */

/*   printf("extended_euclidean algorithm over ints...\n"); */

/*   int d = extended_euclidean_algorithm(a,b,&x,&y); */

/*   printf("a=%i, b=%i, x=%i, y=%i,d=%i\n",a,b,x,y,d); */

/*   mpz_t mpz_a,mpz_b,mpz_x,mpz_y; */
/*   //mpz_init_set_ui(mpz_a,7); */
/*   //mpz_init_set_ui(mpz_b,5); */
/*   mpz_init_set_str(mpz_a,"2199025364993",10); */
/*   mpz_init_set_str(mpz_b,"2199025512449",10); */
/*   mpz_init(mpz_x); */
/*   mpz_init(mpz_y); */

/*   d = mpz_extended_euclidean_algorithm(mpz_a,mpz_b,mpz_x,mpz_y); */
/*   printf("extended_euclidean algorithm over mpzs...\n"); */
  
/*   char *s_a,*s_b,*s_x,*s_y; */
  
/*   s_a = mpz_get_str(NULL,10,mpz_a); */
/*   s_b = mpz_get_str(NULL,10,mpz_b); */
/*   s_x = mpz_get_str(NULL,10,mpz_x); */
/*   s_y = mpz_get_str(NULL,10,mpz_y); */
  
/*   printf("a=%s, b=%s, x=%s, y=%s,d=%i\n",s_a,s_b,s_x,s_y,d); */
/* } */

// given a,b; it returns (d,x,y) such that d = x*a+y*b
// d is returned as value
int extended_euclidean_algorithm(int a,int b,int* ptr_x,int* ptr_y){

  int prevx = 1, x=0,
    prevy=0, y=1;

  while(b!=0){
    int q=a/b;
    int swap;
    swap=x; x=prevx-q*x; prevx=swap;
    swap=y; y=prevy-q*y; prevy=swap;
    swap=a; a=b; b=swap-q*b;
  }

  (*ptr_x)=prevx;
  (*ptr_y)=prevy;

  return a;
  }

int mpz_extended_euclidean_algorithm(mpz_t mpz_a, mpz_t mpz_b, mpz_t ptr_x, mpz_t ptr_y){

  mpz_t prevx, x,
    prevy, y;

  mpz_t a,b;
  mpz_init_set(a,mpz_a);
  mpz_init_set(b,mpz_b);
  
  mpz_init_set_ui(prevx,1);
  mpz_init_set_ui(x,0);
  mpz_init_set_ui(y,1);
  mpz_init_set_ui(prevy,0);

  mpz_t q,r,swap;
  mpz_init(q);
  mpz_init(r);
  mpz_init(swap);
  while(mpz_cmp_ui(b,0)!=0){
    mpz_mdivmod(q,r,a,b);//a/b = q + r/b
    mpz_mul(swap,q,x);
    mpz_sub(prevx,prevx,swap);
    mpz_swap(prevx,x);
    mpz_mul(swap,q,y);
    mpz_sub(prevy,prevy,swap);
    mpz_swap(prevy,y);
    mpz_mul(swap,q,b);
    mpz_sub(a,a,swap);
    mpz_swap(a,b);
  }

  mpz_set(ptr_x,prevx);
  mpz_set(ptr_y,prevy);
  
  return mpz_get_ui(a);//as must be 1, if I get something different -> ERROR
}
