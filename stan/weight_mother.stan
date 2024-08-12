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

array [N, A] real weight;
array [N, A] real height;

array [N] int birthorder;
array [A] real mid_height_by_age;

array [N] int dob; 

array [N] int male;
array [N] int twin;

array [N] int mother_id;
array [N] int father_id;

array [N, A] int father_dead;
array [N, A] int mother_dead;
array [N, A] int mother_unmarried;
array [N, A] int mother_married_to_notfather;
array [N, A] int mother_married_to_father_with_cowife;

array [N, A] int skip;

}

transformed data {

  array [N, A] real log_weight;
  array [N, A] real log_height;

  for (n in 1:N) {
    for (a in 1:A) {

      if (weight[n, a] != -99){
        log_weight[n, a] = log(weight[n, a]);
      }

    }
  }
  
  for (n in 1:N) {
    for (a in 1:A) {

      if (height[n, a] != -99){
        log_height[n, a] = log(height[n, a]);
      }
      
    }
  }
  
}

parameters {

real alpha;
real <lower = 0> scale;

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
vector <lower = 0> [10] a_age_tau;
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

scale ~ exponential(1);

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

      for (a in 1:A) {

        if (weight[n, a] != -99) { // if weight is not unknown
          if (height[n, a] != -99) { // if height is not unknown

            if (skip[n, a] == 0) {

                  log_weight[n, a] ~ normal(
                  alpha + 
                  a_bo[birthorder[n]] +
                  a_year[dob[n] + (a-1)] + 
                  a_mother[mother_id[n]] * mother_sigma +
                  a_father[father_id[n]] * father_sigma +
                  a_age[1, a] + 
                  a_age[2, a] * male[n] + 
                  a_age[3, a] * twin[n] + 
                  a_age[4, a] * father_dead[n, a] +
                  a_age[5, a] * mother_dead[n, a] +
                  a_age[6, a] * mother_unmarried[n, a] +
                  a_age[7, a] * mother_married_to_notfather[n, a] +
                  a_age[8, a] * mother_married_to_father_with_cowife[n, a] +
                  a_age[10, a] * log_height[n, a],
                  scale);

                }

          if (skip[n, a] == 1) {

            log_weight[n, a] ~ normal(
                  alpha + 
                  a_bo[birthorder[n]] +
                  a_year[dob[n] + (a-1)] + 
                  a_mother[mother_id[n]] * mother_sigma +
                  a_father[father_id[n]] * father_sigma +
                  a_age[1, a] + 
                  a_age[2, a] * male[n] + 
                  a_age[3, a] * twin[n] + 
                  a_age[9, a] +
                  a_age[10, a] * log_height[n, a],
                  scale);

            }
          }

        }
      }

  } // n

}

generated quantities {

  array [2] vector [A] m_base;
  array [2] vector [A] m_mother_dead;
  array [2] vector [A] m_mother_unmarried;
  array [2] vector [A] m_mother_married_to_notfather;
  array [2] vector [A] m_mother_married_to_father_with_cowife;
  array [2] vector [A] m_unknown_parent;

  real sum_parent_sigma;
  sum_parent_sigma = mother_sigma + father_sigma;

  for (a in 1:A) {

    for (i in 1:2) { // male or female

      // median predictions for a child born in year id 61, non-twin, other parent alive, of avg height

      m_base[i, a] = exp(alpha +
                         a_bo[1] +
                         a_year[16 + (a-1)] +
                         a_age[1, a] +
                         a_age[2, a] * (i-1) +
                         a_age[10, a] * log(mid_height_by_age[a]));

      m_mother_dead[i, a] = exp(alpha + 
                                a_bo[1] + 
                                a_year[16 + (a-1)] +
                                a_age[1, a] +
                                a_age[2, a] * (i-1) + // male or female offset (when i == 1, female, when i == 2 is male)
                                a_age[5, a] * 1 +
                                a_age[10, a] * log(mid_height_by_age[a]));

      m_mother_unmarried[i, a] = exp(alpha + 
                                     a_bo[1] + 
                                     a_year[16 + (a-1)] +
                                     a_age[1, a] +
                                     a_age[2, a] * (i-1) +
                                     a_age[6, a] * 1 +
                                     a_age[10, a] * log(mid_height_by_age[a]));

      m_mother_married_to_notfather[i, a] = exp(alpha + 
                                                a_bo[1] + 
                                                a_year[16 + (a-1)] +
                                                a_age[1, a] +
                                                a_age[2, a] * (i-1) +
                                                a_age[7, a] * 1 +
                                                a_age[10, a] * log(mid_height_by_age[a]));

      m_mother_married_to_father_with_cowife[i, a] = exp(alpha + 
                                                         a_bo[1] + 
                                                         a_year[16 + (a-1)] +
                                                         a_age[1, a] +
                                                         a_age[2, a] * (i-1) +
                                                         a_age[8, a] * 1 +
                                                         a_age[10, a] * log(mid_height_by_age[a]));

      m_unknown_parent[i, a] = exp(alpha + 
                                   a_bo[1] + 
                                   a_year[16 + (a-1)] +
                                   a_age[1, a] +
                                   a_age[2, a] * (i-1) +
                                   a_age[9, a] * 1 +
                                   a_age[10, a] * log(mid_height_by_age[a]));

    }

  }

}
