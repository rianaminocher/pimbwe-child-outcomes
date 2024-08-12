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

  vector cumulative_product(vector X) {

    return exp(cumulative_sum(log(X)));

  }

}

data {

int N;
int A;
int B;
int Y;
int I;
int L;

array [N, A] int alive;

array [N] int birthorder;
array [N] int dob; 

array [N] int male;
array [N] int twin;

array [N] int mother_id;
array [N] int father_id;

array [N, A] int mother_dead;
array [N, A] int father_dead;
array [N, A] int father_unmarried;
array [N, A] int father_married_to_notmother_monogamy;
array [N, A] int father_married_to_notmother_polygyny;
array [N, A] int father_married_to_mother_polygyny;

array [N, A] int skip;

}

parameters {

real alpha;
real <lower = 0, upper = 1> p_male; // prob missing sex is male
real miss_birth_order; // missing birthorder param

real <lower = 0> mother_sigma;
vector [I] a_mother;

real <lower = 0> father_sigma;
vector [L] a_father;

vector [B] a_bo_raw;
real <lower = 0, upper = 1> a_bo_kappa;
real <lower = 0> a_bo_tau;
real <lower = 0> a_bo_delta;

vector [Y] a_year_raw;
real <lower = 0, upper = 1> a_year_kappa;
real <lower = 0> a_year_tau;
real <lower = 0> a_year_delta;

array [10] vector [A] a_age_raw;
vector <lower = 0, upper = 1> [10] a_age_kappa;
vector <lower = 0> [10] a_age_tau ;
vector <lower = 0> [10] a_age_delta;

}

transformed parameters {

vector [B] a_bo;
vector [Y] a_year;
array [10] vector [A] a_age;

a_bo = GP(B, a_bo_kappa, a_bo_tau, a_bo_delta) * a_bo_raw;
a_bo[16] = miss_birth_order;

a_year = GP(Y, a_year_kappa, a_year_tau, a_year_delta) * a_year_raw;

for (i in 1:10) {
  a_age[i] = GP(A, a_age_kappa[i], a_age_tau[i], a_age_delta[i]) * a_age_raw[i];
}

}

model {

alpha ~ normal(2, 2);
p_male ~ beta(2, 2);

miss_birth_order ~ normal(0, 5);

mother_sigma ~ exponential(1);
a_mother ~ normal(0, 1);

father_sigma ~ exponential(1);
a_father ~ normal(0, 1);

a_bo_raw ~ normal(0, 1);
a_bo_kappa ~ beta(12, 2);
a_bo_tau ~ exponential(1);
a_bo_delta ~ exponential(1);

a_year_raw ~ normal(0, 1);
a_year_kappa ~ beta(12, 2);
a_year_tau ~ exponential(1);
a_year_delta ~ exponential(1);

for (i in 1:10) {

  a_age_raw[i] ~ normal(0, 1);
  a_age_kappa[i] ~ beta(12, 2);
  a_age_tau[i] ~ exponential(1);
  a_age_delta[i] ~ exponential(1);

}

for (n in 1:N) {

  // if sex is not missing
  if (male[n] != -99) {

  male[n] ~ bernoulli(p_male);

  // when a == 1, we cannot condition on survival in (a-1)
    if (alive[n, 1] != -99) {

      if (skip[n, 1] == 0) {

        alive[n, 1] ~ bernoulli_logit(
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * male[n] + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * mother_dead[n, 1] +
        a_age[5, 1] * father_dead[n, 1] +
        a_age[6, 1] * father_unmarried[n, 1] +
        a_age[7, 1] * father_married_to_notmother_monogamy[n, 1] +
        a_age[8, 1] * father_married_to_notmother_polygyny[n, 1] +
        a_age[9, 1] * father_married_to_mother_polygyny[n, 1]);
      }

      if (skip[n, 1] == 1) {

        alive[n, 1] ~ bernoulli_logit(
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * male[n] + 
        a_age[3, 1] * twin[n] + 
        a_age[10, 1]);

      }

    }

      for (a in 2:A) {

        if (alive[n, a] != -99) { // if survival is not unknown
          if (alive[n, a-1] == 1) { // if survived previous year

            if (skip[n, a] == 0) {

                  alive[n, a] ~ bernoulli_logit(
                  alpha +
                  a_bo[birthorder[n]] +
                  a_year[dob[n] + (a-1)] + 
                  a_mother[mother_id[n]] * mother_sigma +
                  a_father[father_id[n]] * father_sigma +
                  a_age[1, a] + 
                  a_age[2, a] * male[n] + 
                  a_age[3, a] * twin[n] + 
                  a_age[4, a] * mother_dead[n, a] +
                  a_age[5, a] * father_dead[n, a] +
                  a_age[6, a] * father_unmarried[n, a] +
                  a_age[7, a] * father_married_to_notmother_monogamy[n, a] +
                  a_age[8, a] * father_married_to_notmother_polygyny[n, a] +
                  a_age[9, a] * father_married_to_mother_polygyny[n, a]);

            }

            if (skip[n, a] == 1) {

              alive[n, a] ~ bernoulli_logit(
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n] + (a-1)] +
              a_mother[mother_id[n]] * mother_sigma +
              a_father[father_id[n]] * father_sigma +
              a_age[1, a] + 
              a_age[2, a] * male[n] + 
              a_age[3, a] * twin[n] + 
              a_age[10, a]);

            }

          }
        }
      }

    } // male known

  // if sex is missing
  // p(y) = p(x==1)p(y|x==1) + p(x==0)p(y|x==0)

  if (male[n] == -99) {

    if (alive[n, 1] != -99) {
      if (skip[n, 1] == 0) {

      real theta_0; // prob given female
      real theta_1; // prob given male

      theta_0 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 0 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * mother_dead[n, 1] +
        a_age[5, 1] * father_dead[n, 1] +
        a_age[6, 1] * father_unmarried[n, 1] +
        a_age[7, 1] * father_married_to_notmother_monogamy[n, 1] +
        a_age[8, 1] * father_married_to_notmother_polygyny[n, 1] +
        a_age[9, 1] * father_married_to_mother_polygyny[n, 1]);

      theta_1 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 1 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * mother_dead[n, 1] +
        a_age[5, 1] * father_dead[n, 1] +
        a_age[6, 1] * father_unmarried[n, 1] +
        a_age[7, 1] * father_married_to_notmother_monogamy[n, 1] +
        a_age[8, 1] * father_married_to_notmother_polygyny[n, 1] +
        a_age[9, 1] * father_married_to_mother_polygyny[n, 1]);

      target += log_mix(p_male, theta_0, theta_1);

    }

      if (skip[n, 1] == 1) {

      real theta_0; // prob given female
      real theta_1; // prob given male

      theta_0 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 0 + 
        a_age[3, 1] * twin[n] + 
        a_age[10, 1]);

      theta_1 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n] + 0] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 1 + 
        a_age[3, 1] * twin[n] + 
        a_age[10, 1]);

      target += log_mix(p_male, theta_0, theta_1);

    }

  }

    for (a in 2:A) {

      if (alive[n, a] != -99) {
        if (alive[n, a-1] == 1) {

          if (skip[n, a] == 0) {

            real theta_0; // prob given female
            real theta_1; // prob given male

            theta_0 = bernoulli_logit_lpmf(
              alive[n, a] |
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n] + (a-1)] + 
              a_mother[mother_id[n]] * mother_sigma +
              a_father[father_id[n]] * father_sigma +
              a_age[1, a] + 
              a_age[2, a] * 0 + 
              a_age[3, a] * twin[n] + 
              a_age[4, a] * mother_dead[n, a] +
              a_age[5, a] * father_dead[n, a] +
              a_age[6, a] * father_unmarried[n, a] +
              a_age[7, a] * father_married_to_notmother_monogamy[n, a] +
              a_age[8, a] * father_married_to_notmother_polygyny[n, a] +
              a_age[9, a] * father_married_to_mother_polygyny[n, a]);

             theta_1 = bernoulli_logit_lpmf(
               alive[n, a] |
               alpha + 
               a_bo[birthorder[n]] + 
               a_year[dob[n] + (a-1)] + 
               a_mother[mother_id[n]] * mother_sigma +
               a_father[father_id[n]] * father_sigma +
               a_age[1, a] + 
               a_age[2, a] * 1 + 
               a_age[3, a] * twin[n] + 
               a_age[4, a] * mother_dead[n, a] +
               a_age[5, a] * father_dead[n, a] +
               a_age[6, a] * father_unmarried[n, a] +
               a_age[7, a] * father_married_to_notmother_monogamy[n, a] +
               a_age[8, a] * father_married_to_notmother_polygyny[n, a] +
               a_age[9, a] * father_married_to_mother_polygyny[n, a]);

            target += log_mix(p_male, theta_1, theta_0);

          }

          if (skip[n, a] == 1) {

            real theta_0; // prob given female
            real theta_1; // prob given male

            theta_0 = bernoulli_logit_lpmf(
              alive[n, a] |
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n] + (a-1)] + 
              a_mother[mother_id[n]] * mother_sigma +
              a_father[father_id[n]] * father_sigma +
              a_age[1, a] + 
              a_age[2, a] * 0 + 
              a_age[3, a] * twin[n] + 
              a_age[10, a]);

             theta_1 = bernoulli_logit_lpmf(
               alive[n, a] |
               alpha + 
               a_bo[birthorder[n]] + 
               a_year[dob[n] + (a-1)] + 
               a_mother[mother_id[n]] * mother_sigma +
               a_father[father_id[n]] * father_sigma +
               a_age[1, a] + 
               a_age[2, a] * 1 + 
               a_age[3, a] * twin[n] + 
               a_age[10, a]);

            target += log_mix(p_male, theta_1, theta_0);

          }

      }
    }

  }

  }//male unknown

  } // n

}

generated quantities {

  array [2] vector [A] p_base;
  array [2] vector [A] p_father_dead;
  array [2] vector [A] p_father_unmarried;
  array [2] vector [A] p_father_married_to_notmother_monogamy;
  array [2] vector [A] p_father_married_to_notmother_polygyny;
  array [2] vector [A] p_father_married_to_mother_polygyny;
  array [2] vector [A] p_unknown_parent;

  array [2] vector [A] m_base;
  array [2] vector [A] m_father_dead;
  array [2] vector [A] m_father_unmarried;
  array [2] vector [A] m_father_married_to_notmother_monogamy;
  array [2] vector [A] m_father_married_to_notmother_polygyny;
  array [2] vector [A] m_father_married_to_mother_polygyny;
  array [2] vector [A] m_unknown_parent;

  real sum_parent_sigma;
  sum_parent_sigma = mother_sigma + father_sigma;

  for (a in 1:A) {

    for (i in 1:2) { // male or female

      p_base[i, a] = inv_logit(alpha +
                               a_bo[1] +
                               a_year[60 + (a-1)] +
                               a_age[1, a] +
                               a_age[2, a] * (i-1));

      p_father_dead[i, a] = inv_logit(alpha + 
                                      a_bo[1] + 
                                      a_year[60 + (a-1)] +
                                      a_age[1, a] +
                                      a_age[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                      a_age[5, a] * 1);

      p_father_unmarried[i, a] = inv_logit(alpha + 
                                           a_bo[1] + 
                                           a_year[60 + (a-1)] +
                                           a_age[1, a] +
                                           a_age[2, a] * (i-1) +
                                           a_age[6, a] * 1);

      p_father_married_to_notmother_monogamy[i, a] = inv_logit(alpha + 
                                                               a_bo[1] + 
                                                               a_year[60 + (a-1)] +
                                                               a_age[1, a] +
                                                               a_age[2, a] * (i-1) +
                                                               a_age[7, a] * 1);

      p_father_married_to_notmother_polygyny[i, a] = inv_logit(alpha + 
                                                               a_bo[1] + 
                                                               a_year[60 + (a-1)] +
                                                               a_age[1, a] +
                                                               a_age[2, a] * (i-1) +
                                                               a_age[8, a] * 1);

      p_father_married_to_mother_polygyny[i, a] = inv_logit(alpha + 
                                                            a_bo[1] + 
                                                            a_year[60 + (a-1)] +
                                                            a_age[1, a] +
                                                            a_age[2, a] * (i-1) +
                                                            a_age[9, a] * 1);

      p_unknown_parent[i, a] = inv_logit(alpha + 
                                         a_bo[1] + 
                                         a_year[60 + (a-1)] +
                                         a_age[1, a] +
                                         a_age[2, a] * (i-1) + 
                                         a_age[10, a]);

    }

  }

  for (i in 1:2) {

    m_base[i] = cumulative_product(p_base[i]);
    m_father_dead[i] = cumulative_product(p_father_dead[i]);
    m_father_unmarried[i] = cumulative_product(p_father_unmarried[i]);
    m_father_married_to_notmother_monogamy[i] = cumulative_product(p_father_married_to_notmother_monogamy[i]);
    m_father_married_to_notmother_polygyny[i] = cumulative_product(p_father_married_to_notmother_polygyny[i]);
    m_father_married_to_mother_polygyny[i] = cumulative_product(p_father_married_to_mother_polygyny[i]);
    m_unknown_parent[i] = cumulative_product(p_unknown_parent[i]);
  }

}
