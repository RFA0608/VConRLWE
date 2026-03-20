#ifndef AUTHENTIC_H
#define AUTHENTIC_H

#include "Group.h"
#include "RLWE.h"

class authentic
{
    private: 
        int poly_degree;
        mpz_class cipher_mod;
        mpz_class g_mod;
        mpz_class g_gen;
        

        // ekf contain(on Group) gamma_0*r_0, alpha_0(r_0*F - r_1), alpha_0*r_0*G, alpha_0*r_0*R,
        //                       gamma_1*r_1, alpha_1(r_1*F - r_0), alpha_1*r_1*G, alpha_1*r_1*R,
        //                       delta*s, beta_0(s*H - r_1), beta_1(s*H - r_0), 
        //                       r_0, r_1, s, respectively. (only memory management)
        std::vector<gvec*> ekf;
        // just pointer mapping (to use easy)
        gvec* g_gamma_0_r_0;
        gvec* g_alpha_0_r_0_F_r_1_m;
        gvec* g_alpha_0_r_0_G;
        gvec* g_alpha_0_r_0_R;
        gvec* g_gamma_1_r_1;
        gvec* g_alpha_1_r_1_F_r_0_m;
        gvec* g_alpha_1_r_1_G;
        gvec* g_alpha_1_r_1_R;
        gvec* g_delta_s;
        gvec* g_beta_0_s_H_r_1_m;
        gvec* g_beta_1_s_H_r_0_m;
        gvec* g_r_0;
        gvec* g_r_1;
        gvec* g_s;

        // key contain r_0, r_1, s, Gamma_0, Gamma_1, alpha_0, alpha_1, beta_0, beta_1, delta, respectively. (only memory management)
        std::vector<poly*> key; 
        // just pointer mapping (to use easy)
        poly* r_0;
        poly* r_1;
        poly* s;
        poly* gamma_0;
        poly* gamma_1;
        poly* alpha_0;
        poly* alpha_1;
        poly* beta_0;
        poly* beta_1;
        poly* delta;

    public:
        authentic(int poly_degree, mpz_class cipher_mod, mpz_class g_mod, mpz_class g_gen):
            poly_degree(poly_degree),
            g_mod(g_mod),
            g_gen(g_gen)
        {
            poly* r_0 = new poly(16 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, r_0);
            this->r_0 = r_0;
            this->key.push_back(r_0);

            poly* r_1 = new poly(16 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, r_1);
            this->r_1 = r_1;
            this->key.push_back(r_1);

            poly* s = new poly(3 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, s);
            this->s = s;
            this->key.push_back(s);

            poly* gamma_0 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, gamma_0);
            this->gamma_0 = gamma_0;
            this->key.push_back(gamma_0);
            
            poly* gamma_1 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, gamma_1);
            this->gamma_1 = gamma_1;
            this->key.push_back(gamma_1);

            poly* alpha_0 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, alpha_0);
            this->alpha_0 = alpha_0;
            this->key.push_back(alpha_0);

            poly* alpha_1 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, alpha_1);
            this->alpha_1 = alpha_1;
            this->key.push_back(alpha_1);

            poly* beta_0 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, beta_0);
            this->beta_0 = beta_0;
            this->key.push_back(beta_0);

            poly* beta_1 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, beta_1);
            this->beta_1 = beta_1;
            this->key.push_back(beta_1);

            poly* delta = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, delta);
            this->delta = delta;
            this->key.push_back(delta);
        };

        ~authentic()
        {
            for(auto& factor : this->ekf)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->ekf.clear();
            for(auto& factor : this->key)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->key.clear();
        }

        std::vector<mr_cipher*> make_encH(cipher* temp, poly* sk, std::vector<std::vector<int64_t>>& P, std::vector<std::vector<int64_t>>& Q)
        {
            std::vector<mr_cipher*> mr_represent(8);
            std::vector<int64_t> pod_matrix(temp->ring_dim, 0LL);
            poly* packed_data = new poly(temp->ring_dim);
            poly* plaintext = new poly(temp->ring_dim);

            for(int i = 0; i < 8; i++)
            {
                if(i < 4)
                {
                    pod_matrix[0] = P[i][0];
                    pod_matrix[1] = P[i][1];

                    batch_encoder::encode(pod_matrix, temp->plain_mod, temp->psi_plain, packed_data);
                    poly_handler::pack_2_plain(packed_data, plaintext);
                    crypto_handler::encrypt(plaintext, sk, temp);
                    
                    mr_represent[i] = new mr_cipher(temp->ring_dim, temp->plain_mod, temp->cipher_mod);
                    format_transform_handler::cipher_2_mr_cipher(temp, mr_represent[i]);
                }
                else
                {
                    pod_matrix[0] = Q[i - 4][0];
                    pod_matrix[1] = 0;

                    batch_encoder::encode(pod_matrix, temp->plain_mod, temp->psi_plain, packed_data);
                    poly_handler::pack_2_plain(packed_data, plaintext);
                    crypto_handler::encrypt(plaintext, sk, temp);

                    mr_represent[i] = new mr_cipher(temp->ring_dim, temp->plain_mod, temp->cipher_mod);
                    format_transform_handler::cipher_2_mr_cipher(temp, mr_represent[i]);
                }
            }

            delete packed_data;
            delete plaintext;

            return mr_represent;
        }

        void make_ekf(std::vector<mr_cipher*> H_enc_mat)
        {
            // ciphertext whole size(N size message part and N size random key part, respectivley)
            int step_size = 2 * this->poly_degree;

            // gamma_0 * r_0 power g
            poly* temp_gamma_0_r_0 = new poly(this->r_0->ring_dim);
            poly_handler::poly_scal_mul_p(this->gamma_0, this->r_0, temp_gamma_0_r_0);
            gvec* gamma_0_r_0 = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_gamma_0_r_0, gamma_0_r_0);
            this->g_gamma_0_r_0 = gamma_0_r_0;
            this->ekf.push_back(gamma_0_r_0);
            delete temp_gamma_0_r_0;

            // alpha_0 * (r_0 * F - r_1) power g
            poly* temp_alpha_0_r_0_F_r_1_m = new poly(this->r_0->ring_dim);
            temp_alpha_0_r_0_F_r_1_m->fill_zero();
            for(int i = 0; i < r_0->ring_dim; i++)
            {
                if(i < step_size)
                {
                    continue;
                }
                else if(i >= step_size && i < 4 * step_size)
                {
                    temp_alpha_0_r_0_F_r_1_m->coeff[i] = this->r_0->coeff[i - step_size];
                }
                else if(i >= 4 * step_size && i < 5 * step_size)
                {
                    continue;
                }
                else
                {
                    temp_alpha_0_r_0_F_r_1_m->coeff[i] = this->r_0->coeff[i - step_size];
                }
            }
            poly* temp_r_1_m1 = new poly(this->r_1->ring_dim);
            poly_handler::poly_neg(this->r_1, temp_r_1_m1);
            poly_handler::poly_add(temp_alpha_0_r_0_F_r_1_m, temp_r_1_m1, temp_alpha_0_r_0_F_r_1_m);
            poly_handler::poly_scal_mul_p(this->alpha_0, temp_alpha_0_r_0_F_r_1_m, temp_alpha_0_r_0_F_r_1_m);
            gvec* alpha_0_r_0_F_r_1_m = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_0_r_0_F_r_1_m, alpha_0_r_0_F_r_1_m);
            this->g_alpha_0_r_0_F_r_1_m = alpha_0_r_0_F_r_1_m;
            this->ekf.push_back(alpha_0_r_0_F_r_1_m);
            delete temp_alpha_0_r_0_F_r_1_m;
            delete temp_r_1_m1;

            // alpha_0 * r_0 * G power g
            poly* temp_alpha_0_r_0_G = new poly(step_size);
            temp_alpha_0_r_0_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_alpha_0_r_0_G->coeff[i] = this->r_0->coeff[i + 3 * step_size];
            }
            poly_handler::poly_scal_mul_p(this->alpha_0, temp_alpha_0_r_0_G, temp_alpha_0_r_0_G);
            gvec* alpha_0_r_0_G = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_0_r_0_G, alpha_0_r_0_G);
            this->g_alpha_0_r_0_G = alpha_0_r_0_G;
            this->ekf.push_back(alpha_0_r_0_G);
            delete temp_alpha_0_r_0_G;

            // alpha_0 * r_0 * R power g
            poly* temp_alpha_0_r_0_R = new poly(step_size);
            temp_alpha_0_r_0_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_alpha_0_r_0_R->coeff[i] = this->r_0->coeff[i + 7 * step_size];
            }
            poly_handler::poly_scal_mul_p(this->alpha_0, temp_alpha_0_r_0_R, temp_alpha_0_r_0_R);
            gvec* alpha_0_r_0_R = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_0_r_0_R, alpha_0_r_0_R);
            this->g_alpha_0_r_0_R = alpha_0_r_0_R;
            this->ekf.push_back(alpha_0_r_0_R);
            delete temp_alpha_0_r_0_R;

             // gamma_1 * r_1 power g
            poly* temp_gamma_1_r_1 = new poly(this->r_1->ring_dim);
            poly_handler::poly_scal_mul_p(this->gamma_1, this->r_1, temp_gamma_1_r_1);
            gvec* gamma_1_r_1 = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_gamma_1_r_1, gamma_1_r_1);
            this->g_gamma_1_r_1 = gamma_1_r_1;
            this->ekf.push_back(gamma_1_r_1);
            delete temp_gamma_1_r_1;

            // alpha_1 * (r_1 * F - r_0) power g
            poly* temp_alpha_1_r_1_F_r_0_m = new poly(this->r_1->ring_dim);
            temp_alpha_1_r_1_F_r_0_m->fill_zero();
            for(int i = 0; i < r_1->ring_dim; i++)
            {
                if(i < step_size)
                {
                    continue;
                }
                else if(i >= step_size && i < 4 * step_size)
                {
                    temp_alpha_1_r_1_F_r_0_m->coeff[i] = this->r_1->coeff[i - step_size];
                }
                else if(i >= 4 * step_size && i < 5 * step_size)
                {
                    continue;
                }
                else
                {
                    temp_alpha_1_r_1_F_r_0_m->coeff[i] = this->r_1->coeff[i - step_size];
                }
            }
            poly* temp_r_0_m1 = new poly(this->r_0->ring_dim);
            poly_handler::poly_neg(this->r_0, temp_r_0_m1);
            poly_handler::poly_add(temp_alpha_1_r_1_F_r_0_m, temp_r_0_m1, temp_alpha_1_r_1_F_r_0_m);
            poly_handler::poly_scal_mul_p(this->alpha_1, temp_alpha_1_r_1_F_r_0_m, temp_alpha_1_r_1_F_r_0_m);
            gvec* alpha_1_r_1_F_r_0_m = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_1_r_1_F_r_0_m, alpha_1_r_1_F_r_0_m);
            this->g_alpha_1_r_1_F_r_0_m = alpha_1_r_1_F_r_0_m;
            this->ekf.push_back(alpha_1_r_1_F_r_0_m);
            delete temp_alpha_1_r_1_F_r_0_m;
            delete temp_r_0_m1;

            // alpha_1 * r_1 * G power g
            poly* temp_alpha_1_r_1_G = new poly(step_size);
            temp_alpha_1_r_1_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_alpha_1_r_1_G->coeff[i] = this->r_1->coeff[i + 3 * step_size];
            }
            poly_handler::poly_scal_mul_p(this->alpha_1, temp_alpha_1_r_1_G, temp_alpha_1_r_1_G);
            gvec* alpha_1_r_1_G = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_1_r_1_G, alpha_1_r_1_G);
            this->g_alpha_1_r_1_G = alpha_1_r_1_G;
            this->ekf.push_back(alpha_1_r_1_G);
            delete temp_alpha_1_r_1_G;

            // alpha_1 * r_1 * R power g
            poly* temp_alpha_1_r_1_R = new poly(step_size);
            temp_alpha_1_r_1_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_alpha_1_r_1_R->coeff[i] = this->r_1->coeff[i + 7 * step_size];
            }
            poly_handler::poly_scal_mul_p(this->alpha_1, temp_alpha_1_r_1_R, temp_alpha_1_r_1_R);
            gvec* alpha_1_r_1_R = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_1_r_1_R, alpha_1_r_1_R);
            this->g_alpha_1_r_1_R = alpha_1_r_1_R;
            this->ekf.push_back(alpha_1_r_1_R);
            delete temp_alpha_1_r_1_R;

            // delta * s power g
            poly* temp_delta_s = new poly(this->s->ring_dim);
            temp_delta_s->fill_zero();
            poly_handler::poly_scal_mul_p(this->delta, this->s, temp_delta_s);
            gvec* delta_s = new gvec(this->s->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_delta_s, delta_s);
            this->g_delta_s = delta_s;
            this->ekf.push_back(delta_s);
            delete temp_delta_s;

            // need to decribe beta_0(s*H - r_1), beta_1(s*H - r_0)
            std::vector<poly*> temp_sH(8);
            for(int i = 0; i < 8; i++)
            {
                temp_sH[i] = new poly(2 * this->poly_degree);
                crypto_handler::pval_mr_mul(this->s, H_enc_mat[i], temp_sH[i]);
            }

            poly* sH = poly_handler::poly_recur_concat(temp_sH);

            poly* temp_r_1_m2 = new poly(this->r_1->ring_dim);
            poly_handler::poly_neg(this->r_1, temp_r_1_m2);
            poly* temp_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_1_m2, temp_s_H_r_1_m);
            poly* temp_beta_0_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_0, temp_s_H_r_1_m, temp_beta_0_s_H_r_1_m);
            gvec* beta_0_s_H_r_1_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_beta_0_s_H_r_1_m, beta_0_s_H_r_1_m);
            this->g_beta_0_s_H_r_1_m = beta_0_s_H_r_1_m;
            this->ekf.push_back(beta_0_s_H_r_1_m);
            delete temp_r_1_m2;
            delete temp_s_H_r_1_m;
            delete temp_beta_0_s_H_r_1_m;
            
            poly* temp_r_0_m2 = new poly(this->r_0->ring_dim);
            poly_handler::poly_neg(this->r_0, temp_r_0_m2);
            poly* temp_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_0_m2, temp_s_H_r_0_m);
            poly* temp_beta_1_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_1, temp_s_H_r_0_m, temp_beta_1_s_H_r_0_m);
            gvec* beta_1_s_H_r_0_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_beta_1_s_H_r_0_m, beta_1_s_H_r_0_m);
            this->g_beta_1_s_H_r_0_m = beta_1_s_H_r_0_m;
            this->ekf.push_back(beta_1_s_H_r_0_m);
            delete temp_r_0_m2;
            delete temp_s_H_r_0_m;
            delete temp_beta_1_s_H_r_0_m;

            for(auto& factor : temp_sH)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            temp_sH.clear();
            delete sH;
            
            // r_0 power g
            gvec* r_0 = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->r_0, r_0);
            this->g_r_0 = r_0;

            // r_1 power g
            gvec* r_1 = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->r_1, r_1);
            this->g_r_1 = r_1;

            // s power g
            gvec* s = new gvec(this->s->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->s, s);
            this->g_s = s;
        }
};

#endif