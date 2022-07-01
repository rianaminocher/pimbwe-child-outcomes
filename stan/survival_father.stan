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

int N;          // total number of individuals
int A;          // number of age categories
int B;          // number of birth order categories
int Y;          // number of calendar years

int alive [N, A];

int birthorder [N];
int dob [N]; 

int male [N];
int twin [N];

int motherdead [N, A]; // whether mother is dead
int fatherdead [N, A]; // whether father is dead
int unmarried [N, A]; // whether father is unmarried
int married_not_to_mother [N, A]; // whether father is married not to mother
int married_to_not_mother [N, A]; // whether father is married to not mother, in addition to mother

int skip [N, A];

}

parameters {

real alpha;
real <lower = 0, upper = 1> p_male; // prob missing sex is male
real miss_birth_order; // missing birthorder param

vector [B] a_bo_raw;
real <lower = 0, upper = 1> a_bo_kappa;
real <lower = 0> a_bo_tau;
real <lower = 0> a_bo_delta;

vector [Y] a_year_raw;
real <lower = 0, upper = 1> a_year_kappa;
real <lower = 0> a_year_tau;
real <lower = 0> a_year_delta;

vector [A] a_age_raw [8];
real <lower = 0, upper = 1> a_age_kappa [8];
real <lower = 0> a_age_tau [8];
real <lower = 0> a_age_delta [8];

}

transformed parameters {

vector [B] a_bo;
vector [Y] a_year;
vector [A] a_age [8];

a_bo = GP(B, a_bo_kappa, a_bo_tau, a_bo_delta) * a_bo_raw;
a_bo[16] = miss_birth_order;

a_year = GP(Y, a_year_kappa, a_year_tau, a_year_delta) * a_year_raw;

for (i in 1:8) {
  a_age[i] = GP(A, a_age_kappa[i], a_age_tau[i], a_age_delta[i]) * a_age_raw[i];
}

}

model {

alpha ~ normal(2, 2);
p_male ~ beta(2, 2);

miss_birth_order ~ normal(0, 5);

a_bo_raw ~ normal(0, 1);
a_bo_kappa ~ beta(12, 2);
a_bo_tau ~ exponential(1);
a_bo_delta ~ exponential(1);

a_year_raw ~ normal(0, 1);
a_year_kappa ~ beta(12, 2);
a_year_tau ~ exponential(1);
a_year_delta ~ exponential(1);

for (i in 1:8) {

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
        a_age[1, 1] + 
        a_age[2, 1] * male[n] + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * motherdead[n, 1] +
        a_age[5, 1] * fatherdead[n, 1] +
        a_age[6, 1] * unmarried[n, 1] +
        a_age[7, 1] * married_not_to_mother[n, 1] +
        a_age[8, 1] * married_to_not_mother[n, 1]);

      }
  }

      for (a in 2:A) {

        if (alive[n, a] != -99) { // if survival is not unknown
          if (alive[n, a-1] == 1) { // if survived previous year
            if (skip[n, a] == 0) {

                  alive[n, a] ~ bernoulli_logit(
                  alpha + 
                  a_bo[birthorder[n]] +
                  a_year[dob[n] + (a - 1)] + 
                  a_age[1, a] + 
                  a_age[2, a] * male[n] + 
                  a_age[3, a] * twin[n] + 
                  a_age[4, a] * motherdead[n, a] +
                  a_age[5, a] * fatherdead[n, a] +
                  a_age[6, a] * unmarried[n, a] +
                  a_age[7, a] * married_not_to_mother[n, a] +
                  a_age[8, a] * married_to_not_mother[n, a]);

            }
          }
        }
      }

    } // male known

  // if sex is missing
  // p(y) = p(x==1)p(y|x==1) + p(x==0)p(y|x==0)

  if (male[n] == -99) {

  real theta_0; // prob given female
  real theta_1; // prob given male

    if (alive[n, 1] != -99) {
      if (skip[n, 1] == 0) {

      theta_0 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_age[1, 1] + 
        a_age[2, 1] * 0 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * motherdead[n, 1] +
        a_age[5, 1] * fatherdead[n, 1] +
        a_age[6, 1] * unmarried[n, 1] +
        a_age[7, 1] * married_not_to_mother[n, 1] +
        a_age[8, 1] * married_to_not_mother[n, 1]);

      theta_1 = bernoulli_logit_lpmf(
        alive[n, 1] |
        alpha + 
        a_bo[birthorder[n]] + 
        a_year[dob[n]] + 
        a_age[1, 1] + 
        a_age[2, 1] * 1 + 
        a_age[3, 1] * twin[n] + 
        a_age[4, 1] * motherdead[n, 1] +
        a_age[5, 1] * fatherdead[n, 1] +
        a_age[6, 1] * unmarried[n, 1] +
        a_age[7, 1] * married_not_to_mother[n, 1] +
        a_age[8, 1] * married_to_not_mother[n, 1]);

      target += log_mix(p_male, theta_1, theta_0);

    }
  }

    for (a in 2:A) {

      if (alive[n, a] != -99) {
        if (alive[n, a-1] == 1) {
          if (skip[n, a] == 0) {

            theta_0 = bernoulli_logit_lpmf(
              alive[n, a] |
              alpha + 
              a_bo[birthorder[n]] + 
              a_year[dob[n] + (a - 1)] + 
              a_age[1, a] + 
              a_age[2, a] * 0 + 
              a_age[3, a] * twin[n] + 
              a_age[4, a] * motherdead[n, a] +
              a_age[5, a] * fatherdead[n, a] +
              a_age[6, a] * unmarried[n, a] +
              a_age[7, a] * married_not_to_mother[n, a] +
              a_age[8, a] * married_to_not_mother[n, a]);

             theta_1 = bernoulli_logit_lpmf(
               alive[n, a] |
               alpha + 
               a_bo[birthorder[n]] + 
               a_year[dob[n] + (a - 1)] + 
               a_age[1, a] + 
               a_age[2, a] * 1 + 
               a_age[3, a] * twin[n] + 
               a_age[4, a] * motherdead[n, a] +
               a_age[5, a] * fatherdead[n, a] +
               a_age[6, a] * unmarried[n, a] +
               a_age[7, a] * married_not_to_mother[n, a] +
               a_age[8, a] * married_to_not_mother[n, a]);

            target += log_mix(p_male, theta_0, theta_1);

          }
        }
      }

    } 
    }//male unknown

  } // n

}

generated quantities {

  vector [A] p_base [2];
  vector [A] p_father_dead [2];
  vector [A] p_father_unmarried [2];
  vector [A] p_father_married_not_to_mother [2];
  vector [A] p_father_married_to_not_mother [2];

  vector [A] cumulative_p_base [2];
  vector [A] cumulative_p_father_dead [2];
  vector [A] cumulative_p_father_unmarried [2];
  vector [A] cumulative_p_father_married_not_to_mother [2];
  vector [A] cumulative_p_father_married_to_not_mother [2];

  for (a in 1:A) {

    for (i in 1:2) { // male or female

      p_base[i, a] = inv_logit(alpha +
                               a_bo[1] +
                               a_year[61 - (a - 1)] +
                               a_age[1, a] +
                               a_age[2, a] * (i-1));

      p_father_dead[i, a] = inv_logit(alpha + 
                                      a_bo[1] + 
                                      a_year[61 - (a - 1)] +
                                      a_age[1, a] +
                                      a_age[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                      a_age[5, a] * 1);

      p_father_unmarried[i, a] = inv_logit(alpha + 
                                           a_bo[1] + 
                                           a_year[61 - (a - 1)] +
                                           a_age[1, a] +
                                           a_age[2, a] * (i-1) +
                                           a_age[6, a] * 1);

      p_father_married_not_to_mother[i, a] = inv_logit(alpha + 
                                                       a_bo[1] + 
                                                       a_year[61 - (a - 1)] +
                                                       a_age[1, a] +
                                                       a_age[2, a] * (i-1) +
                                                       a_age[7, a] * 1);

      p_father_married_to_not_mother[i, a] = inv_logit(alpha + 
                                                       a_bo[1] + 
                                                       a_year[61 - (a - 1)] +
                                                       a_age[1, a] +
                                                       a_age[2, a] * (i-1) +
                                                       a_age[8, a] * 1);
    }

  }

  for (i in 1:2) {

    cumulative_p_base[i] = cumulative_product(p_base[i]);
    cumulative_p_father_dead[i] = cumulative_product(p_father_dead[i]);
    cumulative_p_father_unmarried[i] = cumulative_product(p_father_unmarried[i]);
    cumulative_p_father_married_not_to_mother[i] = cumulative_product(p_father_married_not_to_mother[i]);
    cumulative_p_father_married_to_not_mother[i] = cumulative_product(p_father_married_to_not_mother[i]);

  }


}
