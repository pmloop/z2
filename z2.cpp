// Z2 4D lattice
// cpp adaptation of M. Creutz C code
//

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
# include <random>

using namespace std;

// the lattice is of dimensions SIZE**4
#define SIZE 40

class config {

  private:
  double lattice[SIZE][SIZE][SIZE][SIZE][4];
  public:
  // basic feature of lattice_config
  config ();

  void coldstart();
  void randomstart();
  double staplecal(const int *y, int d);
  double sweep(double beta);

  };

config::config(){
  int x[4];
  int d;
  for (x[0]=0;x[0]<SIZE;x[0]++)
    for (x[1]=0;x[1]<SIZE;x[1]++)
      for (x[2]=0;x[2]<SIZE;x[2]++)
        for (x[3]=0;x[3]<SIZE;x[3]++)
          for (d=0;d<4;d++){
	        lattice[x[0]][x[1]][x[2]][x[3]][d]=1.;}
  return;
}

// utility
// 
// move +1 step along d-th direction
void moveup(int *x, int d) {
  x[d] = x[d] + 1;
  if (x[d]>=SIZE) x[d]=x[d]-SIZE; 
  return;
}

void movedown(int *x, int d) {
  x[d] = x[d] - 1;
  if (x[d]<0) x[d]=x[d]+SIZE; 
  return;
}


void config::coldstart()
  {
  int x[4];
  int d;

  for (x[0]=0;x[0]<SIZE;x[0]++)
    for (x[1]=0;x[1]<SIZE;x[1]++)
      for (x[2]=0;x[2]<SIZE;x[2]++)
        for (x[3]=0;x[3]<SIZE;x[3]++)
          for (d=0;d<4;d++){
	    lattice[x[0]][x[1]][x[2]][x[3]][d]=1.;}
  return;
  }


void config::randomstart()
  {
  int x[4],d;
  // usage: unif(rng) to get a random number
  std::uniform_real_distribution<double> unif(0., 1.);
  std::random_device rand_dev;     // Use random_device to get a random seed.
  std::mt19937_64 rng(rand_dev()); // mt19937 is a good pseudo-random number generator.

  for (x[0]=0; x[0]<SIZE; x[0]++)
    for (x[1]=0; x[1]<SIZE; x[1]++)
      for (x[2]=0; x[2]<SIZE; x[2]++)
        for (x[3]=0; x[3]<SIZE; x[3]++)
          for (d=0; d<4; d++) 
          {
            if ( unif(rng) > 0.5 ){
              lattice[x[0]][x[1]][x[2]][x[3]][d]=1.;
            }
            else{ 
              lattice[x[0]][x[1]][x[2]][x[3]][d]=-1.;
            }
          }
  return;
  }

// staple calculation
double config::staplecal(const int *y, int d)
  {
  double staple,staplesum;    
  int x[4];
  int dperp;

  // copy the address of current site
  x[0] = y[0];
  x[1] = y[1];
  x[2] = y[2];
  x[3] = y[3];

            staplesum=0.;
            for (dperp=0;dperp<4;dperp++){
              if (dperp!=d){
                /*  move around thusly:
                    dperp        6--5
                    ^            |  |
                    |            1--4
                    |            |  |
                    -----> d     2--3  */
                /* plaquette 1234 */
                movedown(x,dperp);
                staple=lattice[x[0]][x[1]][x[2]][x[3]][dperp]
                  *lattice[x[0]][x[1]][x[2]][x[3]][d];
                moveup(x,d);
                staple=staple*lattice[x[0]][x[1]][x[2]][x[3]][dperp];  
                moveup(x,dperp);
                staplesum=staplesum+staple;
                /* plaquette 1456 */
                staple=lattice[x[0]][x[1]][x[2]][x[3]][dperp];
                moveup(x,dperp);
                movedown(x,d);
                staple=staple*lattice[x[0]][x[1]][x[2]][x[3]][d];
                movedown(x,dperp);
                staple=staple*lattice[x[0]][x[1]][x[2]][x[3]][dperp];
                staplesum=staplesum+staple;
              }
	          }

  return staplesum;
  }

  // do one Monte Carlo sweep
  double config::sweep(double beta)
  {

  int x[4],d;
  double staplesum;
  double bplus,bminus;
  double action = 0.0;

  double current_link;

  // usage: unif(rng) to get a random number
  std::uniform_real_distribution<double> unif(0., 1.);
  std::random_device rand_dev;          // Use random_device to get a random seed.
  std::mt19937_64 rng(rand_dev()); // mt19937 is a good pseudo-random number generator.

  for (x[0]=0; x[0]<SIZE; x[0]++)
    for (x[1]=0; x[1]<SIZE; x[1]++)
      for (x[2]=0; x[2]<SIZE; x[2]++)
        for (x[3]=0; x[3]<SIZE; x[3]++)
          for (d=0; d<4; d++) 
          {

            staplesum = staplecal(x, d);

            /* calculate the Boltzmann weight */
            bplus=exp(beta*staplesum);
            bminus=exp(-1.*beta*staplesum);
            bplus=bplus/(bplus+bminus);
            /* the heatbath algorithm */
            if ( unif(rng) < bplus ){
              lattice[x[0]][x[1]][x[2]][x[3]][d]=1.;
	      action +=staplesum;
            }
            else{ 
              lattice[x[0]][x[1]][x[2]][x[3]][d]=-1.;
	      action -=staplesum;
            }
          }
  action /= (SIZE*SIZE*SIZE*SIZE*4*6);
  return 1.-action;
  }


int main(){
  
  // parameters
  double beta_min = 0.;
  double beta_max = 1.;
  const int N_beta = 100;
  double action = 0.;
  double beta, dbeta;

  //static config current;
  static config current;

  dbeta= (beta_max-beta_min)/(N_beta);
  current.coldstart();
  //current.randomstart();
  /* heat it up */
  for (beta=beta_max; beta>beta_min; beta-=dbeta){
    action = current.sweep(beta);
    cout << beta << " " << action  << endl;
  }

  /* cool it down */
  for (beta=beta_min; beta<beta_max; beta+=dbeta){
    action = current.sweep(beta);
    cout << beta << " " << action  << endl;
  }

  return 0;
}

