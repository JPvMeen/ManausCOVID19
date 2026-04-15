functions {
  vector one_strain_corona_model(
      real t, vector y, 
      real beta1, real gamma, real log10_ifr1_intercept, real log10_ifr1_slope,
      real w, real beta1_re, real scalar_ifr1,
      data matrix C,
      data int A, data int J, data int K) {
  
  /** Define the system of ODEs */
  vector[A] ifr1 = pow(10, log10_ifr1_intercept + log10_ifr1_slope * linspaced_vector(A, -A+1, 0));
  vector[A] ifr1_re = scalar_ifr1 * ifr1;
  vector[A] d1 = gamma * ifr1 ./ (1-ifr1);
  vector[A] d1_re = gamma * ifr1_re ./ (1-ifr1_re); 
  
  //States [S, I1, I1_re, R1, W1]
  vector[A] S = y[1 : A];
  vector[J*A] I1_J = y[A+1 : (J+1)*A];
  vector[A] I1_tot = to_matrix(I1_J, A, J) * rep_vector(1.0, J);
  
  vector[J*A] I1_re_J = y[(J+1)*A + 1 : (2*J+1)*A]; // Reinfection stages
  vector[A] I1_re_tot = to_matrix(I1_re_J, A, J) * rep_vector(1.0, J); // Total of I1_re stages
  
  vector[K*A] R1_K = y[(2*J+1)*A + 1 : (2*J+1+K)*A]; // Recovered from Primary infection
  vector[A] R1_tot = to_matrix(R1_K, A, K) * rep_vector(1.0, K); // Total of R1 stages
  vector[A] W1 = y[(2*J+1+K)*A + 1 : (2*J+2+K)*A];  // Waning immunity from R1
  
  //Calculations
  vector[A] N = S + I1_tot + I1_re_tot + R1_tot + W1; 
  vector[A] C_I1_I1_re_N = C * ((I1_tot + I1_re_tot) ./ N); // Added for strain 1
  
  vector[A] incidence1 = C_I1_I1_re_N * beta1 .* S;
  vector[A] incidence1_re = C_I1_I1_re_N * beta1_re .* W1; 
  
  //ODEs as single vector
  vector[A * (2*J + K + 2)] dydt;
  
  dydt[1 : A] = -(incidence1); // dS_dt
  
  // Primary infection (Strain 1)
  dydt[A + 1: 2*A] = incidence1 // first infectious stage 
    - (gamma*J) * I1_J[1: A] 
    - d1 .* I1_J[1: A]; 
  for (j in 2:J) { // other infectious stages
    dydt[j*A + 1: (j+1)*A] = 
    (gamma*J) * I1_J[(j-2)*A + 1: (j-1)*A]
    - (gamma*J) * I1_J[(j-1)*A + 1 : (j)*A] 
    - d1 .* I1_J[(j-1)*A + 1 : (j)*A] ;
  }
  
  // Reinfection (Strain 1)
  dydt[(J+1)*A + 1 : (J+2)*A] = incidence1_re // first I1 reinfectious stage
    - (gamma*J) * I1_re_J[1 : A] 
    - d1_re .* I1_re_J[1 : A]; 
  for (j in 2:J) { // other infectious stages
    dydt[(J+1+j-1)*A + 1: (J+2+j-1)*A] =  
    (gamma*J) * I1_re_J[(j-2)*A + 1 : (j-1)*A]  
    - (gamma*J) * I1_re_J[(j-1)*A + 1: j*A]
    - d1_re .* I1_re_J[(j-1)*A + 1: j*A];
  }
  
  // First R1_K compartment: receives recovery from last infectious stage
  dydt[(2*J+1)*A + 1 : (2*J+2)*A] =
    (gamma*J) .* (I1_J[(J-1)*A + 1 : J*A] + I1_re_J[(J-1)*A + 1 : J*A])
    - w*K .* R1_K[1:A];

  // Remaining R1_K compartments
  for (k in 2:K) {
    dydt[(2*J+1+k-1)*A + 1 : (2*J+2+k-1)*A] =
      w*K .* R1_K[(k-2)*A + 1: (k-1)*A]
      - w*K .* R1_K[(k-1)*A + 1: k*A];
  }
    
  dydt[(2*J+1+K)*A + 1 : (2*J+2+K)*A] = w*K * R1_K[(K-1)*A + 1 : K*A] - incidence1_re; // dW1_dt
  
  return dydt;
}

vector two_strain_corona_model(
      real t, vector y, 
      real beta1, real gamma, real log10_ifr1_intercept, real log10_ifr1_slope,
      real w, real beta1_re, real scalar_ifr1,
      real beta2, real log10_ifr2_intercept, real log10_ifr2_slope,
      real beta2_re, real scalar_ifr2,
      data matrix C,
      data int A, data int J, data int K) {
  
  /** Define the system of ODEs */

  vector[A] ifr1 = pow(10, log10_ifr1_intercept + log10_ifr1_slope * linspaced_vector(A, -A+1, 0));
  vector[A] ifr1_re = scalar_ifr1 * ifr1;
  vector[A] d1 = gamma * ifr1 ./ (1-ifr1);
  vector[A] d1_re = gamma * ifr1_re ./ (1-ifr1_re); 
  
  vector[A] ifr2 = pow(10, log10_ifr2_intercept + log10_ifr2_slope * linspaced_vector(A, -A+1, 0));
  vector[A] ifr2_re = scalar_ifr2 * ifr2;
  vector[A] d2 = gamma * ifr2 ./ (1-ifr2);
  vector[A] d2_re = gamma * ifr2_re ./ (1-ifr2_re); // Death rate for reinfection (strain 2)
  
  //States [S, I1, I1_re, R1, W1, I2, I2_re, R2, W2]
  vector[A] S = y[1 : A];
  vector[J*A] I1_J = y[A+1 : (J+1)*A];
  vector[A] I1_tot = to_matrix(I1_J, A, J) * rep_vector(1.0, J);
  
  vector[J*A] I1_re_J = y[(J+1)*A + 1 : (2*J+1)*A]; // Reinfection stages
  vector[A] I1_re_tot = to_matrix(I1_re_J, A, J) * rep_vector(1.0, J); // Total of I1_re stages

  vector[K*A] R1_K = y[(2*J+1)*A + 1 : (2*J+1+K)*A]; // Recovered from Primary infection
  vector[A] R1_tot = to_matrix(R1_K, A, K) * rep_vector(1.0, K); // Total of R1 stages
  vector[A] W1 = y[(2*J+1+K)*A + 1 : (2*J+2+K)*A];  // Waning immunity from R1
  
  vector[J*A] I2_J = y[(2*J+2+K)*A + 1 : (3*J+2+K)*A];
  vector[J*A] I2_re_J = y[(3*J+2+K)*A + 1 : (4*J+2+K)*A];

  vector[A] I2_tot = to_matrix(I2_J, A, J) * rep_vector(1.0, J);
  vector[A] I2_re_tot = to_matrix(I2_re_J, A, J) * rep_vector(1.0, J);
  
  vector[K*A] R2_K = y[(4*J+2+K)*A + 1 : (4*J+2+2*K)*A];
  vector[A] R2_tot = to_matrix(R2_K, A, K) * rep_vector(1.0, K);
  vector[A] W2 = y[(4*J+2+2*K)*A + 1 : (4*J+3+2*K)*A];
  
  //Calculations
  vector[A] N = S + I1_tot + I1_re_tot + R1_tot + W1 + I2_tot + I2_re_tot + R2_tot + W2; 
  vector[A] C_I1_I1_re_N = C * ((I1_tot + I1_re_tot) ./ N); 
  vector[A] C_I2_I2_re_N = C * ((I2_tot + I2_re_tot) ./ N);
  
  vector[A] incidence1 = C_I1_I1_re_N * beta1 .* S;
  vector[A] incidence1_re = C_I1_I1_re_N * beta1_re .* W1; 
  
  vector[A] incidence2 = C_I2_I2_re_N * beta2 .* S; 
  vector[A] incidence2_re = C_I2_I2_re_N * beta2_re .* W1; 
  
  //ODEs as single vector
  vector[(4*J + 2*K + 3)*A] dydt;
  
  dydt[1 : A] = -incidence1 -incidence2; // dS_dt
  
  // Primary infection (Strain 1)
  dydt[A + 1: 2*A] = incidence1 // first infectious stage 
    - (gamma*J) * I1_J[1: A] 
    - d1 .* I1_J[1: A]; 
  for (j in 2:J) { // other infectious stages
    dydt[j*A + 1: (j+1)*A] = 
    (gamma*J) * I1_J[(j-2)*A + 1: (j-1)*A]
    - (gamma*J) * I1_J[(j-1)*A + 1 : (j)*A] 
    - d1 .* I1_J[(j-1)*A + 1 : (j)*A] ;
  }
  
  // Reinfection (Strain 1)
  dydt[(J+1)*A + 1 : (J+2)*A] = incidence1_re // first I1 reinfectious stage
    - (gamma*J) * I1_re_J[1 : A] 
    - d1_re .* I1_re_J[1 : A]; 
  for (j in 2:J) { // other infectious stages
    dydt[(J+1+j-1)*A + 1: (J+2+j-1)*A] =  
    (gamma*J) * I1_re_J[(j-2)*A + 1 : (j-1)*A]  
    - (gamma*J) * I1_re_J[(j-1)*A + 1: j*A]
    - d1_re .* I1_re_J[(j-1)*A + 1: j*A];
  }
  
  // First R1_K compartment: receives recovery from last infectious stage
  dydt[(2*J+1)*A + 1 : (2*J+2)*A] =
    (gamma*J) .* (I1_J[(J-1)*A + 1 : J*A] + I1_re_J[(J-1)*A + 1 : J*A])
    - w*K .* R1_K[1:A];

  // Remaining R1_K compartments
  for (k in 2:K) {
    dydt[(2*J+1+k-1)*A + 1 : (2*J+2+k-1)*A] =
      w*K .* R1_K[(k-2)*A + 1: (k-1)*A]
      - w*K .* R1_K[(k-1)*A + 1: k*A];
  }
  
  dydt[(2*J+1+K)*A + 1 : (2*J+2+K)*A] = w*K * R1_K[(K-1)*A + 1 : K*A] - incidence1_re - incidence2_re; // dW1_dt
  
  // Primary infection (Strain 2)
  dydt[(2*J+2+K)*A + 1 : (2*J+3+K)*A] = incidence2 // first I2 infectious stage
    - (gamma*J) * I2_J[1 : A] 
    - d2 .* I2_J[1 : A]; 
  for (j in 2:J) { // other infectious stages
    dydt[(2*J+2+K + j-1)*A + 1 : (2*J+3+K + j-1)*A] =  
      (gamma*J) * I2_J[(j-2)*A + 1 : (j-1)*A]  
      - (gamma*J) * I2_J[(j-1)*A + 1: j*A]
      - d2 .* I2_J[(j-1)*A + 1: j*A];
  }
  
  // Reinfection (Strain 2)
  dydt[(3*J+2+K)*A + 1 : (3*J+3+K)*A] = incidence2_re // first I2 reinfectious stage
    - (gamma*J) * I2_re_J[1 : A]
    - d2_re .* I2_re_J[1 : A]; 
  for (j in 2:J) { // other infectious stages
    dydt[(3*J+2+K + j-1)*A + 1 : (3*J+3+K + j-1)*A] =  
      (gamma*J) * I2_re_J[(j-2)*A + 1 : (j-1)*A]  
      - (gamma*J) * I2_re_J[(j-1)*A + 1: j*A]
      - d2_re .* I2_re_J[(j-1)*A + 1: j*A];
  }
  
  // First R2_K compartment: receives recovery from I2 and I2_re
  dydt[(4*J+2+K)*A + 1 : (4*J+3+K)*A] =
    (gamma*J) .* (I2_J[(J-1)*A + 1 : J*A] + I2_re_J[(J-1)*A + 1 : J*A])
    - w*K .* R2_K[1:A];

  // Remaining R2_K compartments
  for (k in 2:K) {
    dydt[(4*J + 2 + K + k-1)*A + 1 : (4*J + 3 + K + k-1)*A] =
      w*K .* R2_K[(k-2)*A + 1 : (k-1)*A]
      - w*K .* R2_K[(k-1)*A + 1 : k*A];
  }
  
  dydt[(4*J+2+2*K)*A + 1 : (4*J+3+2*K)*A] =
    w*K * R2_K[(K-1)*A + 1 : K*A]; //dW2_dt
  
  return dydt;
}

array[] vector solve_odes(
  real t0, vector y0, array[] real ts, 
  data real rel_tol, data real abs_tol, int max_num_steps,
  real t1, real beta1, real gamma, 
  real log10_ifr1_intercept, real log10_ifr1_slope, 
  real w, real beta1_re, real scalar_ifr1, 
  real beta2, real log10_ifr2_intercept, real log10_ifr2_slope,
  real beta2_re, real scalar_ifr2, 
  real inoculum_2, vector demography,
  data matrix Cunp,
  data int A, data int J, data int K) {
  
  int num_ts = num_elements(ts);
  int r = rank(append_row(t1+1, to_vector(ts)), 1); // Integrate first wave up to t1
  
  int dim1 = A*(2*J + K + 2);      // S, I1_J, I1_re_J, R1_K, W1
  int dim  = A*(4*J + 2*K + 3);    // full system

  array[num_ts] vector[dim] y = rep_array(rep_vector(0.0, dim), num_ts);
  
  // Integrate from 0 to t1 using one-strain model
  array[r] vector[dim1] ys1 = ode_adams_tol(one_strain_corona_model, 
    y0, t0, ts[1:r], 
    rel_tol, abs_tol, max_num_steps,
    beta1, gamma, log10_ifr1_intercept, log10_ifr1_slope,
    w, beta1_re, scalar_ifr1, 
    Cunp,
    A, J, K); 
  
  y[1:r, 1:dim1] = ys1; // Fill in y up to t1
  
  // Create starting state for strain 2
  vector[dim] y1 = y[r]; // Take final output of ys1
  y1[1:A] -= (inoculum_2/2)*(demography/sum(demography)); // Adjust S
  y1[(2*J + K + 1)*A + 1 : (2*J + K + 2)*A] -= (inoculum_2/2)*(demography/sum(demography)); // Adjust W1
  y1[(2*J + K + 2)*A + 1: (2*J + K + 2 + J)*A] = to_vector(rep_matrix((inoculum_2/2)*(demography/sum(demography))/J, J)); // I2
  y1[(2*J + K + 2 + J)*A + 1 : (2*J + K + 2 + 2*J)*A] = to_vector(rep_matrix((inoculum_2/2)*(demography/sum(demography))/J, J)); // I2_re
  
  // Integrate second wave using two-strain model
  array[num_ts-r] vector[dim] ys2 = ode_adams_tol(two_strain_corona_model, 
    y1, t1, ts[r+1:num_ts], 
    rel_tol, abs_tol, max_num_steps,
    beta1, gamma, log10_ifr1_intercept, log10_ifr1_slope,
    w, beta1_re, scalar_ifr1, 
    beta2, log10_ifr2_intercept, log10_ifr2_slope,
    beta2_re, scalar_ifr2,
    Cunp,
    A, J, K); 
    
  y[r+1:num_ts, ] = ys2;
  
  return y; // Return combined results
}

}


data {
  /* preliminaries */
  int<lower=1> A; // Number of Ageclasses
  int<lower=1> J; // number of compartments for Erlang-distributed infectious period 
  int<lower=1> K; // number of compartments for Erlang-distributed recovered period 
  real<upper=0> t0; //start integrator here with introduction of strain 1
  real<lower=0> t1; //introduction time strain 2

  /* contact matrices */
  matrix[A, A] Cunp; // unperturbed contact matrix

  /* demography */
  vector[A] demography; // demographic composition of Manaus 2020

  /* integration times */ 
  int<lower=1> numdays_all;
  array[numdays_all] int<lower=1> ts_all; // should include at least all the timepoints in ts_hospdeaths and ts_sero

  /* hospitaldeath data of SIVEP-GRIPE by death date  */
  int<lower=1> numdays_hospdeaths; //number of excess death datapoints
  array[numdays_hospdeaths] int ts_hospdeaths; //time datapoints
  array[numdays_hospdeaths, A] int hospdeaths;  //excess death datapoints at each time datapoint 

  /* serological data from Manaus */
  int<lower=1> numdays_sero;
  array[numdays_sero] int ts_sero; //time datapoints
  array[numdays_sero, (A-2)] int sero_num_sampled;  //sample datapoints at each time datapoint 
  array[numdays_sero, (A-2)] int sero_num_pos;  //positive sample datapoints at each time datapoint 
  array[numdays_sero, (A-2)] int t_sero_start;  //for each ts_sero what what the survey start for that agegroup
  array[numdays_sero, (A-2)] int t_sero_end;  //for each ts_sero what what the survey end for that agegroup

  /* serological data from Manaus */
  int<lower=1> numdays_phylo;
  array[numdays_phylo] int ts_phylo; //time datapoints
  array[numdays_phylo] int phylo_num_sampled;  //sample datapoints at each time datapoint 
  array[numdays_phylo] int phylo_num_pos;  //positive sample datapoints at each time datapoint 

  /* ODE integrator settings */
  real<lower = 0> rel_tol;
  real<lower = 0> abs_tol;
  int max_num_steps;
}

transformed data {
  real<lower=-3, upper=0> log10_ifr1_intercept = -1.1828; //For 70+ the mean age is ~78 (based on midpoints), which means that the IFR should be 6.5%
}

parameters {
  real<lower=log(1), upper=log(1e4)> log_inoculum_1; // initial infections B.1 at time t0, between 1-10k
  real<lower=log(1), upper=log(1e4)> log_inoculum_2; // initial infections P.1 at time t1, between 1-10k
  real<lower=0, upper=1> beta1; //  infection rate B.1
  real<lower=0, upper=1> beta2; // infection rate P.1
  real<lower=0, upper=1> scalar_beta1_re; // waned immunity reduced susceptiblity against B.1 scalar
  real<lower=0, upper=1> delta_betare;  //Additional term for reduced susceptiblity against P.1 

  real<lower=0, upper=1> gamma; //recovery rate

  real<lower=0> log10_ifr1_slope;

  real<lower=-3,upper=0> log10_ifr2_intercept; 
  real<lower=0, upper=1> log10_ifr2_slope;

  // Raw parameters
  real<lower=0, upper=100> scalar_ifr1_raw; 
  real<lower=0, upper=100> delta_ifr;
  real<lower=0, upper=100> w_raw; 

  // overdispersion
  real<lower=0> inv_phi_hosp;

}

transformed parameters {
  real inoculum_1 = exp(log_inoculum_1);
  real inoculum_2 = exp(log_inoculum_2);

  vector[A * (2*J + K + 2)] y0; // Includes [S, I1, I1_re, R1, W1]

  y0[1:A] = demography - inoculum_1 * (demography / sum(demography)); // S
  y0[A + 1: (J+1)*A] = to_vector(rep_matrix(inoculum_1 * (demography / sum(demography)) / J, J)); // I1
  y0[(J+1)*A + 1: (2*J+1)*A] = rep_vector(0.0, A * J); // I1_re (Reinfection for strain 1)
  y0[(2*J+1)*A + 1: (2*J+1+K)*A] = rep_vector(0.0, A * K); // R1
  y0[(2*J+1+K)*A + 1: (2*J+2+K)*A] = rep_vector(0.0, A); // W1

  // Scaled parameters
  real<lower=0, upper=1> beta1_re = scalar_beta1_re * beta1; // Reinfection rate B.1

  real<lower=0, upper=1> scalar_beta2_re = scalar_beta1_re + (1 - scalar_beta1_re) * delta_betare; // waned immunity reduced susceptiblity against P.1 scalar
  real<lower=0, upper=1> beta2_re = scalar_beta2_re * beta2; // Reinfection rate P.1

  real<lower=0, upper=1> scalar_ifr1 = scalar_ifr1_raw * 0.01; // Reinfection fatality scalar B.1
  real<lower=0, upper=100> scalar_ifr2_raw = scalar_ifr1_raw + (100-scalar_ifr1_raw) * (delta_ifr*0.01); // Raw reinfection fatality scalar P.1
  real<lower=0, upper=1> scalar_ifr2 = scalar_ifr2_raw * 0.01; // Reinfection fatality scalar P.1 
  real<lower=0, upper=1> w = w_raw * 0.01; // Immunity waning rate

  /* IFR by age class */
  vector<lower=0>[A] ifr1 = pow(10, log10_ifr1_intercept + log10_ifr1_slope * linspaced_vector(A, -A+1, 0)); // Age-specific ifr B.1 
  vector<lower=0>[A] ifr1_re = scalar_ifr1 * ifr1; // Age-specific reinfection fatality rate B.1 
  vector<lower=0>[A] ifr2 = pow(10, log10_ifr2_intercept + log10_ifr2_slope * linspaced_vector(A, -A+1, 0)); // Age-specific ifr P.1 
  vector<lower=0>[A] ifr2_re = scalar_ifr2 * ifr2; // Age-specific reinfection fatality rate P.1 

  // Death rates by age class
  vector[A] d1 = gamma * (ifr1 ./ (1 - ifr1)); // Age-specific death rate B.1 
  vector[A] d1_re = gamma * (ifr1_re ./ (1 - ifr1_re)); // Age-specific reinfection death rate B.1 
  vector[A] d2 = gamma * (ifr2 ./ (1 - ifr2)); // Age-specific death rate P.1 
  vector[A] d2_re = gamma * (ifr2_re ./ (1 - ifr2_re)); // Age-specific reinfection death rate P.1 

  // overdispersion
  real<lower=0> phi_hosp = 1.0 / inv_phi_hosp; 

  /* For easy reference */
  matrix[numdays_all, A] Susceptible;
  matrix[numdays_all, A] Infected1 = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Infected1_re = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Recovered1 = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Waned1;
  matrix[numdays_all, A] Infected2 = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Infected2_re = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Recovered2 = rep_matrix(0.0, numdays_all, A);
  matrix[numdays_all, A] Waned2;
  matrix[numdays_all, A] total_N;

  array[numdays_all] vector[(4*J + 2*K + 3) * A] y_hat;  // full dim = A*(4J + 2K + 3)

  /* integrate ODEs and take states */  
  y_hat = solve_odes(
      t0, y0, ts_all, 
      rel_tol, abs_tol, max_num_steps,
      t1, beta1, gamma, 
      log10_ifr1_intercept, log10_ifr1_slope, 
      w, beta1_re, scalar_ifr1,
      beta2, log10_ifr2_intercept, log10_ifr2_slope,
      beta2_re, scalar_ifr2, 
      inoculum_2, demography,
      Cunp,
      A, J, K);

  /* extract trajectories */
  for (i in 1:numdays_all) {
    Susceptible[i, :] = transpose(y_hat[i][1:A]); // from 1 to A
      
    for (j in 1:J) {
     Infected1[i, :] += transpose(y_hat[i][j*A + 1:(j+1)*A]); // I1
     Infected1_re[i, :] += transpose(y_hat[i][(J+j)*A + 1:(J+j+1)*A]); // I1_re
    }

    Recovered1[i, :] = rep_row_vector(0.0, A);
    for (k in 1:K) {
     Recovered1[i, :] += transpose(y_hat[i][(2*J + k)*A + 1 : (2*J + k + 1)*A]); // R1_k
    }

    Waned1[i, :] = transpose(y_hat[i][(2*J + K + 1)*A + 1 : (2*J + K + 2)*A]); // W1

    for (j in 1:J) {
      Infected2[i, :] += transpose(y_hat[i][(2*J + K + 1 + j)*A + 1 : (2*J + K + 1 + j + 1)*A]); // I2
      Infected2_re[i, :] += transpose(y_hat[i][(2*J + K + 1 + J + j)*A + 1 : (2*J + K + 1 + J + j + 1)*A]); // I2_re
    }

    Recovered2[i, :] = rep_row_vector(0.0, A);
    for (k in 1:K) {
      Recovered2[i, :] += transpose(y_hat[i][(2*J + K + 1 + 2*J + k)*A + 1 : (2*J + K + 1 + 2*J + k + 1)*A]); // R2_k
    }
    Waned2[i, :] = transpose(y_hat[i][(2*J + K + 2 + 2*J + K)*A + 1 : (2*J + K + 2 + 2*J + K + 1)*A]); // W2

    total_N[i, :] = Susceptible[i, :] 
                + Infected1[i, :] + Infected1_re[i, :]
                + Infected2[i, :] + Infected2_re[i, :]
                + Recovered1[i, :] + Recovered2[i, :]
                + Waned1[i, :] + Waned2[i, :];
}


  /* likelihood contributions */
  matrix[numdays_hospdeaths, A] log_likes_hospdeaths; // log-likelihood contributions of hospitalisations  
  matrix[numdays_sero, A - 2] log_likes_sero;  // log-likelihood contributions of serological data
  vector[numdays_phylo] log_likes_phylo;

  /* expected hospdeaths */
  matrix[numdays_hospdeaths, A] expected_hospdeaths;   // for credible intervals

  for (i in 1:numdays_hospdeaths) {
      int idx = ts_hospdeaths[i];
      vector[A] lambda1 = d1 .* transpose(Infected1[idx, :]);
      vector[A] lambda1_re = d1_re .* transpose(Infected1_re[idx, :]); // Added reinfection deaths
      vector[A] lambda2 = d2 .* transpose(Infected2[idx, :]);
      vector[A] lambda2_re = d2_re .* transpose(Infected2_re[idx, :]);
      vector[A] lambda = lambda1 + lambda1_re + lambda2 + lambda2_re; // Include I1_re deaths
      expected_hospdeaths[i, :] = transpose(lambda);

      for (a in 1:A) {
          log_likes_hospdeaths[i,a] = neg_binomial_2_lpmf(hospdeaths[i, a] | fmax(expected_hospdeaths[i, a], 1e-10), phi_hosp);
      }
  }

  /* sero */
  matrix[numdays_sero, A - 2] p_sero_pos;

  for (i in 1:numdays_sero) {
      for(a in 1:(A-2)) {
        int survey_start = t_sero_start[i,a];
        int survey_end = t_sero_end[i,a];
        int survey_length = survey_end - survey_start + 1;
        real accum = 0.0;

        for(t in survey_start:survey_end) {
          accum += (Recovered1[t, a+1] + Recovered2[t, a+1]) / total_N[t, a+1];
        }

        p_sero_pos[i, a] = accum / (survey_length);
        log_likes_sero[i, a] = binomial_lpmf(sero_num_pos[i, a] | sero_num_sampled[i, a], fmax(p_sero_pos[i, a], 1e-10)); 
      }
  }

  /* phylo */
  vector[numdays_phylo] p_phylo_pos;

  for (i in 1:numdays_phylo) {
      int idx = ts_phylo[i];
      p_phylo_pos[i] = sum(Infected2[idx, ] + Infected2_re[idx, ]) / 
        sum(Infected1[idx, ] + Infected1_re[idx, ] + Infected2[idx, ] + Infected2_re[idx, ]);

      log_likes_phylo[i] = binomial_lpmf(phylo_num_pos[i] | phylo_num_sampled[i], fmax(p_phylo_pos[i], 1e-10)); 
      
  }

}

model {
  //priors
  scalar_beta1_re ~ beta(2,1); // (95% prior coverage: 0.16-0.99)
  delta_betare ~ beta(1, 2); // vague prior, we assume that the reinfection protection with Gamma P1 is similar to non-P.1

  gamma ~ inv_gamma(22.6, 2.44); // 1/infectious period; 95% prior coverage 4.2-15 days

  //log10_ifr1_intercept ~ normal(-1.1828, 0.07); //For 70+ the mean age is ~78 (based on midpoints), which   means that the IFR should be 6.5%
  log10_ifr1_slope ~ normal(0.524, 0.013); // slope 0.524 * age, we use 10 year agebins
  
  log10_ifr2_intercept ~ normal(-1.1828, 0.3); //we assume that for the oldest agegroup the ifr should be similar

  // raw parameters with jacobian corrections
  w_raw ~ inv_gamma(2.40, 0.0153/0.01); // 1/waning period; 95% prior coverage 25 - 409 days, E(T) = 157 days
  target += beta_lpdf(scalar_ifr1_raw / 100 | 1, 2) - log(100); //vague prior, we assume reinfected individuals rarely die (95% prior coverage: 0.01-0.84)
  target += beta_lpdf(delta_ifr / 100 | 1, 2) - log(100); // vague prior, we assume that the reinfection protection with Gamma P1 is similar to non-P.1 (95% prior coverage: 0.01-0.84)

  // Overdispersion
  inv_phi_hosp ~ lognormal(-6, 3); // vague prior

  // Likelihood contributions
  target += sum(log_likes_hospdeaths);
  target += sum(log_likes_sero);
  target += sum(log_likes_phylo);
}

generated quantities {
// Log-likelihood calculations
real log_lik = sum(log_likes_hospdeaths) + sum(log_likes_sero) + sum(log_likes_phylo);
array[A*numdays_hospdeaths + (A - 2)*numdays_sero + 1*numdays_phylo] real log_lik_vec = 
append_array(
   append_array(to_array_1d(log_likes_hospdeaths), to_array_1d(log_likes_sero)),
   to_array_1d(log_likes_phylo)
  );

/* Additional Epidemiological Quantities */
real R0_1 = beta1 / gamma;
real R0_2 = beta2 / gamma;
real recovery_time = 1 / gamma;
real waning_time = 1 / w;
  
/* Simulated Hospital Deaths */
matrix<lower=0>[numdays_hospdeaths, A] simulated_hospdeaths;  // For prediction intervals
  
for (i in 1:numdays_hospdeaths) {  
    simulated_hospdeaths[i,:] = to_row_vector(neg_binomial_2_rng(expected_hospdeaths[i,:], phi_hosp));
}

/* Serological Data */
matrix<lower=0>[numdays_sero, A-2] expected_serodata;   // For credible intervals
matrix<lower=0>[numdays_sero, A-2] simulated_serodata;  // For prediction intervals 

for (i in 1:numdays_sero) {
    expected_serodata[i,:] = to_row_vector(sero_num_sampled[i,:]) .* p_sero_pos[i,:];
    simulated_serodata[i,:] = to_row_vector(binomial_rng(sero_num_sampled[i,:], p_sero_pos[i,:]));
}


/* Genomic Data */
vector<lower=0>[numdays_phylo] expected_phylodata;   // For credible intervals
vector<lower=0>[numdays_phylo] simulated_phylodata;  // For prediction intervals 

for (i in 1:numdays_phylo) {
  expected_phylodata[i] = phylo_num_sampled[i] * p_phylo_pos[i]; 
  simulated_phylodata[i] = binomial_rng(phylo_num_sampled[i], p_phylo_pos[i]);
}

}
