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

int alive [N, A];

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

vector [A] a_age_raw [9];
real <lower = 0, upper = 1> a_age_kappa [9];
real <lower = 0> a_age_tau [9];
real <lower = 0> a_age_delta [9];

}

transformed parameters {

vector [B] a_bo;
vector [Y] a_year;
vector [A] a_age [9];

a_bo = GP(B, a_bo_kappa, a_bo_tau, a_bo_delta) * a_bo_raw;
a_bo[16] = miss_birth_order;

a_year = GP(Y, a_year_kappa, a_year_tau, a_year_delta) * a_year_raw;

for (i in 1:9) {
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

for (i in 1:9) {

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
        a_year[dob[n]] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * male[n] + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * father_dead[n, 1] +
        a_age[5, 1] * mother_dead[n, 1] +
        a_age[6, 1] * mother_unmarried[n, 1] +
        a_age[7, 1] * mother_married_to_notfather[n, 1] +
        a_age[8, 1] * mother_married_to_father_with_cowife[n, 1]);

      }

      if (skip[n, 1] == 1) {

        alive[n, 1] ~ bernoulli_logit(
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_age[1, 1] + 
        a_age[2, 1] * male[n] + 
        a_age[3, 1] * twin[n] + 
        a_age[9, 1]);

      }

    }

      for (a in 2:A) {

        if (alive[n, a] != -99) { // if survival is not unknown
          if (alive[n, a-1] == 1) { // if survived previous year

            if (skip[n, a] == 0) {

                  alive[n, a] ~ bernoulli_logit(
                  alpha + 
                  a_bo[birthorder[n]] +
                  a_year[dob[n]] + 
                  a_mother[mother_id[n]] * mother_sigma +
                  a_father[father_id[n]] * father_sigma +
                  a_age[1, a] + 
                  a_age[2, a] * male[n] + 
                  a_age[3, a] * twin[n] + 
                  a_age[4, a] * father_dead[n, a] +
                  a_age[5, a] * mother_dead[n, a] +
                  a_age[6, a] * mother_unmarried[n, a] +
                  a_age[7, a] * mother_married_to_notfather[n, a] +
                  a_age[8, a] * mother_married_to_father_with_cowife[n, a]);
                }

            if (skip[n, 1] == 1) {

              alive[n, a] ~ bernoulli_logit(
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n]] + 
              a_age[1, a] + 
              a_age[2, a] * male[n] + 
              a_age[3, a] * twin[n] + 
              a_age[9, a]);

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
        a_year[dob[n]] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 0 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * father_dead[n, 1] +
        a_age[5, 1] * mother_dead[n, 1] +
        a_age[6, 1] * mother_unmarried[n, 1] +
        a_age[7, 1] * mother_married_to_notfather[n, 1] + 
        a_age[8, 1] * mother_married_to_father_with_cowife[n, 1]);

      theta_1 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_mother[mother_id[n]] * mother_sigma +
        a_father[father_id[n]] * father_sigma +
        a_age[1, 1] + 
        a_age[2, 1] * 1 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * father_dead[n, 1] +
        a_age[5, 1] * mother_dead[n, 1] +
        a_age[6, 1] * mother_unmarried[n, 1] +
        a_age[7, 1] * mother_married_to_notfather[n, 1] + 
        a_age[8, 1] * mother_married_to_father_with_cowife[n, 1]);

      target += log_mix(p_male, theta_1, theta_0);

  }

      if (skip[n, 1] == 1) {

      real theta_0; // prob given female
      real theta_1; // prob given male

      theta_0 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_age[1, 1] + 
        a_age[2, 1] * 0 + 
        a_age[3, 1] * twin[n] + 
        a_age[9, 1]);

      theta_1 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_age[1, 1] + 
        a_age[2, 1] * 1 + 
        a_age[3, 1] * twin[n] + 
        a_age[9, 1]);

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
              a_year[dob[n]] + 
              a_mother[mother_id[n]] * mother_sigma +
              a_father[father_id[n]] * father_sigma +
              a_age[1, a] + 
              a_age[2, a] * 0 + 
              a_age[3, a] * twin[n] + 
              a_age[4, a] * father_dead[n, a] +
              a_age[5, a] * mother_dead[n, a] +
              a_age[6, a] * mother_unmarried[n, a] +
              a_age[7, a] * mother_married_to_notfather[n, a] +
              a_age[8, a] * mother_married_to_father_with_cowife[n, a]);

             theta_1 = bernoulli_logit_lpmf(
               alive[n, a] |
               alpha + 
               a_bo[birthorder[n]] + 
               a_year[dob[n]] + 
               a_mother[mother_id[n]] * mother_sigma +
               a_father[father_id[n]] * father_sigma +
               a_age[1, a] + 
               a_age[2, a] * 1 + 
               a_age[3, a] * twin[n] + 
               a_age[4, a] * father_dead[n, a] +
               a_age[5, a] * mother_dead[n, a] +
               a_age[6, a] * mother_unmarried[n, a] +
               a_age[7, a] * mother_married_to_notfather[n, a] +
               a_age[8, a] * mother_married_to_father_with_cowife[n, a]);

              target += log_mix(p_male, theta_0, theta_1);

        }

          if (skip[n, a] == 1) {

            real theta_0; // prob given female
            real theta_1; // prob given male

            theta_0 = bernoulli_logit_lpmf(
              alive[n, a] |
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n]] + 
              a_age[1, a] + 
              a_age[2, a] * 0 + 
              a_age[3, a] * twin[n] + 
              a_age[9, a]);

             theta_1 = bernoulli_logit_lpmf(
               alive[n, a] |
               alpha + 
               a_bo[birthorder[n]] + 
               a_year[dob[n]] + 
               a_age[1, a] + 
               a_age[2, a] * 1 + 
               a_age[3, a] * twin[n] + 
               a_age[9, a]);

            target += log_mix(p_male, theta_1, theta_0);

          }

      }
    }

    } 

    }//male unknown

  } // n

}

generated quantities {

  vector [A] p_base [2];
  vector [A] p_mother_dead [2];
  vector [A] p_mother_unmarried [2];
  vector [A] p_mother_married_to_notfather [2];
  vector [A] p_mother_married_to_father_with_cowife [2];
  vector [A] p_unknown_parent [2];

  vector [A] m_base [2];
  vector [A] m_mother_dead [2];
  vector [A] m_mother_unmarried [2];
  vector [A] m_mother_married_to_notfather [2];
  vector [A] m_mother_married_to_father_with_cowife [2];
  vector [A] m_unknown_parent [2];

  real sum_parent_sigma;
  sum_parent_sigma = mother_sigma + father_sigma;

  for (a in 1:A) {

    for (i in 1:2) { // male or female

      p_base[i, a] = inv_logit(alpha +
                               a_bo[1] +
                               a_year[60] +
                               a_age[1, a] +
                               a_age[2, a] * (i-1));

      p_mother_dead[i, a] = inv_logit(alpha + 
                                      a_bo[1] + 
                                      a_year[60] +
                                      a_age[1, a] +
                                      a_age[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                      a_age[5, a] * 1);

      p_mother_unmarried[i, a] = inv_logit(alpha + 
                                           a_bo[1] + 
                                           a_year[60] +
                                           a_age[1, a] +
                                           a_age[2, a] * (i-1) +
                                           a_age[6, a] * 1);

      p_mother_married_to_notfather[i, a] = inv_logit(alpha + 
                                                      a_bo[1] + 
                                                      a_year[60] +
                                                      a_age[1, a] +
                                                      a_age[2, a] * (i-1) +
                                                      a_age[7, a] * 1);

      p_mother_married_to_father_with_cowife[i, a] = inv_logit(alpha + 
                                                               a_bo[1] + 
                                                               a_year[60] +
                                                               a_age[1, a] +
                                                               a_age[2, a] * (i-1) +
                                                               a_age[8, a] * 1);

      p_unknown_parent[i, a] = inv_logit(alpha + 
                                         a_bo[1] + 
                                         a_year[60] +
                                         a_age[1, a] +
                                         a_age[2, a] * (i-1) + 
                                         a_age[9, a]);

  }

  for (i in 1:2) {

    m_base[i] = cumulative_product(p_base[i]);
    m_mother_dead[i] = cumulative_product(p_mother_dead[i]);
    m_mother_unmarried[i] = cumulative_product(p_mother_unmarried[i]);
    m_mother_married_to_notfather[i] = cumulative_product(p_mother_married_to_notfather[i]);
    m_mother_married_to_father_with_cowife[i] = cumulative_product(p_mother_married_to_father_with_cowife[i]);
    m_unknown_parent[i] = cumulative_product(p_unknown_parent[i]);

  }

}

}
