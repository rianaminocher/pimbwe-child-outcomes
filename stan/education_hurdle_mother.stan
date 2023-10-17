functions {

  matrix GP(int A, real kappa, real tau, real delta) {

    matrix[A, A] omega;

    for (a in 1:(A-1)) {
    for (b in (a+1):A) {
    omega[a, b] = kappa * exp(-tau * ((b-a)^2 / A^2));
    omega[b, a] = omega[a, b];
    }
    }

    for (i in 1:A) {
    omega[i, i] = 1;
    }

    return delta * cholesky_decompose(omega);
  }

}

data {

int N;
int A;
int B;
int Y;
int I;
int L;

int edu [N, A];

int birthorder [N];
int dob [N]; 

int male [N];
int twin [N];

int mother_id [N];
int father_id [N];

int father_dead [N, A];
int mother_dead [N, A];
int mother_unmarried [N, A];
int mother_married_to_notfather [N, A];
int mother_married_to_father_with_cowife [N, A];

int skip [N, A];

}

parameters {

real alpha_t;
real alpha_l;

real miss_birth_order; // missing birthorder param

real <lower = 0> mother_sigma_t;
real <lower = 0> mother_sigma_l;

vector [I] a_mother_t;
vector [I] a_mother_l;

real <lower = 0> father_sigma_t;
real <lower = 0> father_sigma_l;

vector [L] a_father_t;
vector [L] a_father_l;

vector [B] a_bo_raw_t;
vector [B] a_bo_raw_l;

real <lower = 0, upper = 1> a_bo_kappa_t;
real <lower = 0, upper = 1> a_bo_kappa_l;
real <lower = 0> a_bo_tau_t;
real <lower = 0> a_bo_tau_l;
real <lower = 0> a_bo_delta_t;
real <lower = 0> a_bo_delta_l;

vector [Y] a_year_raw_t;
vector [Y] a_year_raw_l;

real <lower = 0, upper = 1> a_year_kappa_t;
real <lower = 0, upper = 1> a_year_kappa_l;
real <lower = 0> a_year_tau_t;
real <lower = 0> a_year_tau_l;
real <lower = 0> a_year_delta_t;
real <lower = 0> a_year_delta_l;

vector [A] a_age_raw_t [9];
vector [A] a_age_raw_l [9];

real <lower = 0, upper = 1> a_age_kappa_t [9];
real <lower = 0, upper = 1> a_age_kappa_l [9];
real <lower = 0> a_age_tau_t [9];
real <lower = 0> a_age_tau_l [9];
real <lower = 0> a_age_delta_t [9];
real <lower = 0> a_age_delta_l [9];

}

transformed parameters {

vector [B] a_bo_t;
vector [B] a_bo_l;

vector [Y] a_year_t;
vector [Y] a_year_l;

vector [A] a_age_t [9];
vector [A] a_age_l [9];

a_bo_t = GP(B, a_bo_kappa_t, a_bo_tau_t, a_bo_delta_t) * a_bo_raw_t;
a_bo_l = GP(B, a_bo_kappa_l, a_bo_tau_l, a_bo_delta_l) * a_bo_raw_l;

a_bo_t[16] = miss_birth_order;
a_bo_l[16] = miss_birth_order;

a_year_t = GP(Y, a_year_kappa_t, a_year_tau_t, a_year_delta_t) * a_year_raw_t;
a_year_l = GP(Y, a_year_kappa_l, a_year_tau_l, a_year_delta_l) * a_year_raw_l;

for (i in 1:9) {
  a_age_t[i] = GP(A, a_age_kappa_t[i], a_age_tau_t[i], a_age_delta_t[i]) * a_age_raw_t[i];
  a_age_l[i] = GP(A, a_age_kappa_l[i], a_age_tau_l[i], a_age_delta_l[i]) * a_age_raw_l[i];
}

}

model {

alpha_t ~ normal(2, 2);
alpha_l ~ normal(2, 2);

miss_birth_order ~ normal(0, 5);

mother_sigma_t ~ exponential(1);
a_mother_t ~ normal(0, 1);

mother_sigma_l ~ exponential(1);
a_mother_l ~ normal(0, 1);

father_sigma_t ~ exponential(1);
a_father_t ~ normal(0, 1);

father_sigma_l ~ exponential(1);
a_father_l ~ normal(0, 1);

a_bo_raw_t ~ normal(0, 1);
a_bo_kappa_t ~ beta(12, 2);
a_bo_tau_t ~ exponential(1);
a_bo_delta_t ~ exponential(1);

a_bo_raw_l ~ normal(0, 1);
a_bo_kappa_l ~ beta(12, 2);
a_bo_tau_l ~ exponential(1);
a_bo_delta_l ~ exponential(1);

a_year_raw_t ~ normal(0, 1);
a_year_kappa_t ~ beta(12, 2);
a_year_tau_t ~ exponential(1);
a_year_delta_t ~ exponential(1);

a_year_raw_l ~ normal(0, 1);
a_year_kappa_l ~ beta(12, 2);
a_year_tau_l ~ exponential(1);
a_year_delta_l ~ exponential(1);

for (i in 1:9) {

  a_age_raw_t[i] ~ normal(0, 1);
  a_age_kappa_t[i] ~ beta(12, 2);
  a_age_tau_t[i] ~ exponential(1);
  a_age_delta_t[i] ~ exponential(1);

  a_age_raw_l[i] ~ normal(0, 1);
  a_age_kappa_l[i] ~ beta(12, 2);
  a_age_tau_l[i] ~ exponential(1);
  a_age_delta_l[i] ~ exponential(1);

}

for (n in 1:N) {

      for (a in 1:A) {

        if (edu[n, a] != -99) {

            if (skip[n, a] == 0) {

            	real theta;
            	real lambda;

            	theta = inv_logit(
            	alpha_t + 
            	a_bo_t[birthorder[n]] +
            	a_year_t[dob[n] + (a-1)] +
            	a_mother_t[mother_id[n]] * mother_sigma_t +
            	a_father_t[father_id[n]] * father_sigma_t +
            	a_age_t[1, a] + 
            	a_age_t[2, a] * male[n] + 
            	a_age_t[3, a] * twin[n] + 
            	a_age_t[4, a] * father_dead[n, a] +
            	a_age_t[5, a] * mother_dead[n, a] + 
            	a_age_t[6, a] * mother_unmarried[n, a] +
            	a_age_t[7, a] * mother_married_to_notfather[n, a] +
            	a_age_t[8, a] * mother_married_to_father_with_cowife[n, a]);

            	lambda = exp(
            	alpha_l + 
            	a_bo_l[birthorder[n]] +
            	a_year_l[dob[n] + (a-1)] +
            	a_mother_l[mother_id[n]] * mother_sigma_l +
            	a_father_l[father_id[n]] * father_sigma_l +
            	a_age_l[1, a] + 
            	a_age_l[2, a] * male[n] + 
            	a_age_l[3, a] * twin[n] + 
              a_age_l[4, a] * father_dead[n, a] +
              a_age_l[5, a] * mother_dead[n, a] + 
              a_age_l[6, a] * mother_unmarried[n, a] +
              a_age_l[7, a] * mother_married_to_notfather[n, a] +
              a_age_l[8, a] * mother_married_to_father_with_cowife[n, a]);

            	if (edu[n, a] == 0) {

              target += log_sum_exp(log(theta), log1m(theta) + poisson_lpmf(edu[n, a] | lambda)); 

              } else {

            	    target += log1m(theta) + poisson_lpmf(edu[n, a] | lambda);
            	}

            } // skip

            if (skip[n, a] == 1) {

            		real theta;
            		real lambda;

                theta = inv_logit(
                    alpha_t + 
                    a_bo_t[birthorder[n]] +
                    a_year_t[dob[n] + (a-1)] +
                    a_mother_t[mother_id[n]] * mother_sigma_t +
                    a_father_t[father_id[n]] * father_sigma_t +
                    a_age_t[1, a] + 
                    a_age_t[2, a] * male[n] + 
                    a_age_t[3, a] * twin[n] + 
                    a_age_t[9, a]);

            		lambda = exp(
                 	  alpha_l + 
                  	a_bo_l[birthorder[n]] +
                  	a_year_l[dob[n] + (a-1)] +
                    a_mother_l[mother_id[n]] * mother_sigma_l +
                    a_father_l[father_id[n]] * father_sigma_l +
                  	a_age_l[1, a] + 
                  	a_age_l[2, a] * male[n] + 
                  	a_age_l[3, a] * twin[n] + 
                  	a_age_l[9, a]);

                  	if (edu[n, a] == 0) {

                  	target += log_sum_exp(log(theta), log1m(theta) + poisson_lpmf(edu[n, a] | lambda)); 

                    } else {
 
                  		target += log1m(theta) + poisson_lpmf(edu[n, a] | lambda);

                  	}

            	}

            }

        } // a

  } // n

}

generated quantities {

  vector [A] m_base_t [2];
  vector [A] m_base_l [2];
  vector [A] m_base [2];

  vector [A] m_mother_dead_t [2];
  vector [A] m_mother_dead_l [2];
  vector [A] m_mother_dead [2];

  vector [A] m_mother_unmarried_t [2];
  vector [A] m_mother_unmarried_l [2];
  vector [A] m_mother_unmarried [2];

  vector [A] m_mother_married_to_notfather_t [2];
  vector [A] m_mother_married_to_notfather_l [2];
  vector [A] m_mother_married_to_notfather [2];

  vector [A] m_mother_married_to_father_with_cowife_t [2];
  vector [A] m_mother_married_to_father_with_cowife_l [2];
  vector [A] m_mother_married_to_father_with_cowife [2];

  vector [A] m_unknown_parent_t [2];
  vector [A] m_unknown_parent_l [2];
  vector [A] m_unknown_parent [2];

  real sum_parent_sigma_t;
  sum_parent_sigma_t = mother_sigma_t + father_sigma_t;

  real sum_parent_sigma_l;
  sum_parent_sigma_l = mother_sigma_t + father_sigma_l;

  for (a in 1:A) {

    for (i in 1:2) { // male or female

      m_base_t[i, a] = inv_logit(alpha_t +
                                 a_bo_t[1] +
                                 a_year_t[16 + (a-1)] +
                                 a_age_t[1, a] +
                                 a_age_t[2, a] * (i-1));

      m_base_l[i, a] = exp(alpha_l +
                           a_bo_l[1] +
                           a_year_l[16 + (a-1)] +
                           a_age_l[1, a] +
                           a_age_l[2, a] * (i-1));

      m_base[i, a] = (1 - bernoulli_rng(m_base_t[i, a])) * poisson_rng(m_base_l[i, a]);
 
      m_mother_dead_t[i, a] = inv_logit(alpha_t + 
                                 	      a_bo_t[1] + 
                                        a_year_t[16 + (a-1)] +
                                        a_age_t[1, a] +
                                        a_age_t[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                        a_age_t[5, a] * 1);

      m_mother_dead_l[i, a] = exp(alpha_l + 
                                  a_bo_l[1] + 
                                  a_year_l[16 + (a-1)] +
                                  a_age_l[1, a] +
                                  a_age_l[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                  a_age_l[5, a] * 1);

      m_mother_dead[i, a] = (1 - bernoulli_rng(m_mother_dead_t[i, a])) * poisson_rng(m_mother_dead_l[i, a]);

      m_mother_unmarried_t[i, a] = inv_logit(alpha_t + 
                                     		     a_bo_t[1] + 
                                             a_year_t[16 + (a-1)] +
                                     		     a_age_t[1, a] +
                                     		     a_age_t[2, a] * (i-1) +
                                     		     a_age_t[6, a] * 1);

      m_mother_unmarried_l[i, a] = exp(alpha_l + 
                                       a_bo_l[1] + 
                                       a_year_l[16 + (a-1)] +
                                       a_age_l[1, a] +
                                       a_age_l[2, a] * (i-1) +
                                       a_age_l[6, a] * 1);

      m_mother_unmarried[i, a] = (1 - bernoulli_rng(m_mother_unmarried_t[i, a])) * poisson_rng(m_mother_unmarried_l[i, a]);                                       

      m_mother_married_to_notfather_t[i, a] = inv_logit(alpha_t + 
                                                        a_bo_t[1] + 
                                                        a_year_t[16 + (a-1)] +
                                                        a_age_t[1, a] +
                                                        a_age_t[2, a] * (i-1) +
                                                        a_age_t[7, a] * 1);

      m_mother_married_to_notfather_l[i, a] = exp(alpha_l + 
                                                  a_bo_l[1] + 
                                                  a_year_l[16 + (a-1)] +
                                                  a_age_l[1, a] +
                                                  a_age_l[2, a] * (i-1) +
                                                  a_age_l[7, a] * 1);

      m_mother_married_to_notfather[i, a] = (1 - bernoulli_rng(m_mother_married_to_notfather_t[i, a])) * poisson_rng(m_mother_married_to_notfather_l[i, a]);   

      m_mother_married_to_father_with_cowife_t[i, a] = inv_logit(alpha_t + 
                                                                 a_bo_t[1] + 
                                                                 a_year_t[16 + (a-1)] +
                                                                 a_age_t[1, a] +
                                                                 a_age_t[2, a] * (i-1) +
                                                                 a_age_t[8, a] * 1);

      m_mother_married_to_father_with_cowife_l[i, a] = exp(alpha_l + 
                                                           a_bo_l[1] + 
                                                           a_year_l[16 + (a-1)] +
                                                           a_age_l[1, a] +
                                                           a_age_l[2, a] * (i-1) +
                                                           a_age_l[8, a] * 1);

      m_mother_married_to_father_with_cowife[i, a] = (1 - bernoulli_rng(m_mother_married_to_father_with_cowife_t[i, a])) * poisson_rng(m_mother_married_to_father_with_cowife_l[i, a]);  

      m_unknown_parent_t[i, a] = inv_logit(alpha_t + 
      									                 a_bo_t[1] + 
                                   		   a_year_t[16 + (a-1)] +
                                   		   a_age_t[1, a] +
                                   		   a_age_t[2, a] * (i-1) +
                                   		   a_age_t[9, a]);

      m_unknown_parent_l[i, a] = exp(alpha_l + 
                                   a_bo_l[1] + 
                                   a_year_l[16 + (a-1)] +
                                   a_age_l[1, a] +
                                   a_age_l[2, a] * (i-1) +
                                   a_age_l[9, a]);

      m_unknown_parent[i, a] = (1 - bernoulli_rng(m_unknown_parent_t[i, a])) * poisson_rng(m_unknown_parent_l[i, a]);  

    }

  }

}
