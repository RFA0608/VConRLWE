#ifndef AUTHENTIC_H
#define AUTHENTIC_H

#include "Group.h"
#include "ECC.h"
#include "RLWE.h"

class authentic
{
    public:
        int poly_degree;
        mpz_class cipher_mod;
        mpz_class g_mod;
        mpz_class g_gen;
        

        // ekf contain(on Group) gamma_0*r_0, alpha_0(r_0*F - r_1), (r_0*F - r_1), r_0*G, r_0*R,
        //                       gamma_1*r_1, alpha_1(r_1*F - r_0), (r_1*F - r_0), r_1*G, r_1*R,
        //                       beta_0(s*H - r_1), (s*H - r_1), beta_1(s*H - r_0), (s*H - r_0),
        //                       r_0, r_1, s, respectively. (only memory management)
        std::vector<gvec*> ekf;
        // just pointer mapping (to use easy)
        gvec* g_gamma_0_r_0;
        gvec* g_alpha_0_r_0_F_r_1_m;
        gvec* g_r_0_F_r_1_m;
        gvec* g_r_0_G;
        gvec* g_r_0_R;
        gvec* g_gamma_1_r_1;
        gvec* g_alpha_1_r_1_F_r_0_m;
        gvec* g_r_1_F_r_0_m;
        gvec* g_r_1_G;
        gvec* g_r_1_R;
        gvec* g_beta_0_s_H_r_1_m;
        gvec* g_s_H_r_1_m;
        gvec* g_beta_1_s_H_r_0_m;
        gvec* g_s_H_r_0_m;
        gvec* g_r_0;
        gvec* g_r_1;
        gvec* g_s;

        // key contain r_0, r_1, s, Gamma_0, Gamma_1, alpha_0, alpha_1, beta_0, beta_1, respectively. (only memory management)
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


        std::vector<poly*> mem_stack_new;
        std::vector<poly*> mem_stack_pre;

        poly* arx_state_even;
        poly* arx_state_odd;

        bool timing = false;
        std::vector<mpz_class> even_nu;
        std::vector<mpz_class> even_mu;
        std::vector<mpz_class> odd_nu;
        std::vector<mpz_class> odd_mu;

        poly* Y;
        poly* U_re;
        poly* U;


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

            this->mem_stack_new.resize(16);
            this->mem_stack_pre.resize(16);

            poly* temp;
            for(int i = 0; i < 16; i++)
            {
                temp = new poly(this->poly_degree);
                this->mem_stack_new[i] = temp;
                temp = new poly(this->poly_degree);
                this->mem_stack_pre[i] = temp;
            }

            this->arx_state_even = new poly(this->poly_degree * this->mem_stack_new.size());
            this->arx_state_odd = new poly(this->poly_degree * this->mem_stack_new.size());

            this->even_nu.resize(4);
            this->even_mu.resize(2);
            this->odd_nu.resize(4);
            this->odd_mu.resize(2);

            this->Y = new poly(this->poly_degree * 2);
            this->U_re = new poly(this->poly_degree * 2);
            this->U = new poly(this->poly_degree * 3);
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

            for(auto& factor : this->mem_stack_new)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_stack_new.clear();
            for(auto& factor : this->mem_stack_pre)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_stack_pre.clear();

            delete this->arx_state_even;
            delete this->arx_state_odd;

            this->even_nu.clear();
            this->even_mu.clear();
            this->odd_nu.clear();
            this->odd_mu.clear();

            delete this->Y;
            delete this->U_re;
            delete this->U;
        }
        
        void make_ekf(std::vector<cipher*>& P_enc, std::vector<cipher*>& Q_enc)
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

            // alpha_0 * (r_0 * F - r_1) power g and (r_0 * F - r_1) power g
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
            gvec* r_0_F_r_1_m = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_0_r_0_F_r_1_m, r_0_F_r_1_m);
            this->g_r_0_F_r_1_m = r_0_F_r_1_m;
            poly_handler::poly_scal_mul_p(this->alpha_0, temp_alpha_0_r_0_F_r_1_m, temp_alpha_0_r_0_F_r_1_m);
            gvec* alpha_0_r_0_F_r_1_m = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_0_r_0_F_r_1_m, alpha_0_r_0_F_r_1_m);
            this->g_alpha_0_r_0_F_r_1_m = alpha_0_r_0_F_r_1_m;
            this->ekf.push_back(alpha_0_r_0_F_r_1_m);
            this->ekf.push_back(r_0_F_r_1_m);
            delete temp_alpha_0_r_0_F_r_1_m;
            delete temp_r_1_m1;

            // r_0 * G power g
            poly* temp_r_0_G = new poly(step_size);
            temp_r_0_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_0_G->coeff[i] = this->r_0->coeff[i + 3 * step_size];
            }
            gvec* r_0_G = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_r_0_G, r_0_G);
            this->g_r_0_G = r_0_G;
            this->ekf.push_back(r_0_G);
            delete temp_r_0_G;

            // r_0 * R power g
            poly* temp_r_0_R = new poly(step_size);
            temp_r_0_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_0_R->coeff[i] = this->r_0->coeff[i + 7 * step_size];
            }
            gvec* r_0_R = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_r_0_R, r_0_R);
            this->g_r_0_R = r_0_R;
            this->ekf.push_back(r_0_R);
            delete temp_r_0_R;

            // gamma_1 * r_1 power g
            poly* temp_gamma_1_r_1 = new poly(this->r_1->ring_dim);
            poly_handler::poly_scal_mul_p(this->gamma_1, this->r_1, temp_gamma_1_r_1);
            gvec* gamma_1_r_1 = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_gamma_1_r_1, gamma_1_r_1);
            this->g_gamma_1_r_1 = gamma_1_r_1;
            this->ekf.push_back(gamma_1_r_1);
            delete temp_gamma_1_r_1;

            // alpha_1 * (r_1 * F - r_0) power g and (r_0 * F - r_1) power g
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
            gvec* r_1_F_r_0_m = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_1_r_1_F_r_0_m, r_1_F_r_0_m);
            this->g_r_1_F_r_0_m = r_1_F_r_0_m;
            poly_handler::poly_scal_mul_p(this->alpha_1, temp_alpha_1_r_1_F_r_0_m, temp_alpha_1_r_1_F_r_0_m);
            gvec* alpha_1_r_1_F_r_0_m = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_alpha_1_r_1_F_r_0_m, alpha_1_r_1_F_r_0_m);
            this->g_alpha_1_r_1_F_r_0_m = alpha_1_r_1_F_r_0_m;
            this->ekf.push_back(alpha_1_r_1_F_r_0_m);
            this->ekf.push_back(r_1_F_r_0_m);
            delete temp_alpha_1_r_1_F_r_0_m;
            delete temp_r_0_m1;

            // r_1 * G power g
            poly* temp_r_1_G = new poly(step_size);
            temp_r_1_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_1_G->coeff[i] = this->r_1->coeff[i + 3 * step_size];
            }
            gvec* r_1_G = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_r_1_G, r_1_G);
            this->g_r_1_G = r_1_G;
            this->ekf.push_back(r_1_G);
            delete temp_r_1_G;

            // r_1 * R power g
            poly* temp_r_1_R = new poly(step_size);
            temp_r_1_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_1_R->coeff[i] = this->r_1->coeff[i + 7 * step_size];
            }
            gvec* r_1_R = new gvec(step_size, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_r_1_R, r_1_R);
            this->g_r_1_R = r_1_R;
            this->ekf.push_back(r_1_R);
            delete temp_r_1_R;

            // beta_0(s*H - r_1), (s*H - r_1), beta_1(s*H - r_0), (s*H - r_0), respectively
            // mr_cipher* comp_H_enc_mat =  new mr_cipher(P_enc[0]->ring_dim, P_enc[0]->plain_mod, P_enc[0]->cipher_mod);
            std::vector<poly*> temp_sH(8);
            for(int i = 0; i < 8; i++)
            {
                temp_sH[i] = new poly(2 * this->poly_degree);
                if(i < 4)
                {
                    // format_transform_handler::cipher_2_mr_cipher(P_enc[i], comp_H_enc_mat);
                    // crypto_handler::pval_mr_mul(this->s, comp_H_enc_mat, temp_sH[i]);
                    crypto_handler::pval_mrlike_mul(this->s, P_enc[i], temp_sH[i]);
                }
                else
                {
                    // format_transform_handler::cipher_2_mr_cipher(Q_enc[i - 4], comp_H_enc_mat);
                    // crypto_handler::pval_mr_mul(this->s, comp_H_enc_mat, temp_sH[i]);
                    crypto_handler::pval_mrlike_mul(this->s, Q_enc[i - 4], temp_sH[i]);
                }
            }

            // delete comp_H_enc_mat;

            poly* sH = poly_handler::poly_recur_concat(temp_sH);

            for(auto& factor : temp_sH)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            temp_sH.clear();

            poly* temp_r_1_m2 = new poly(this->r_1->ring_dim);
            poly_handler::poly_neg(this->r_1, temp_r_1_m2);
            poly* temp_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_1_m2, temp_s_H_r_1_m);
            gvec* s_H_r_1_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_s_H_r_1_m, s_H_r_1_m);
            this->g_s_H_r_1_m = s_H_r_1_m;
            poly* temp_beta_0_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_0, temp_s_H_r_1_m, temp_beta_0_s_H_r_1_m);
            gvec* beta_0_s_H_r_1_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_beta_0_s_H_r_1_m, beta_0_s_H_r_1_m);
            this->g_beta_0_s_H_r_1_m = beta_0_s_H_r_1_m;
            this->ekf.push_back(beta_0_s_H_r_1_m);
            this->ekf.push_back(s_H_r_1_m);
            delete temp_r_1_m2;
            delete temp_s_H_r_1_m;
            delete temp_beta_0_s_H_r_1_m;
            
            poly* temp_r_0_m2 = new poly(this->r_0->ring_dim);
            poly_handler::poly_neg(this->r_0, temp_r_0_m2);
            poly* temp_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_0_m2, temp_s_H_r_0_m);
            gvec* s_H_r_0_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_s_H_r_0_m, s_H_r_0_m);
            this->g_s_H_r_0_m = s_H_r_0_m;
            poly* temp_beta_1_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_1, temp_s_H_r_0_m, temp_beta_1_s_H_r_0_m);
            gvec* beta_1_s_H_r_0_m = new gvec(sH->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(temp_beta_1_s_H_r_0_m, beta_1_s_H_r_0_m);
            this->g_beta_1_s_H_r_0_m = beta_1_s_H_r_0_m;
            this->ekf.push_back(beta_1_s_H_r_0_m);
            this->ekf.push_back(s_H_r_0_m);
            delete temp_r_0_m2;
            delete temp_s_H_r_0_m;
            delete temp_beta_1_s_H_r_0_m;
            delete sH;
            
            // r_0 power g
            gvec* r_0 = new gvec(this->r_0->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->r_0, r_0);
            this->g_r_0 = r_0;
            this->ekf.push_back(r_0);

            // r_1 power g
            gvec* r_1 = new gvec(this->r_1->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->r_1, r_1);
            this->g_r_1 = r_1;
            this->ekf.push_back(r_1);

            // s power g
            gvec* s = new gvec(this->s->ring_dim, this->g_mod, this->g_gen);
            group_handler::poly_2_gvec(this->s, s);
            this->g_s = s;
            this->ekf.push_back(s);
        }

        void generate_proof(std::vector<cipher*>& mem_y_new, std::vector<cipher*>& mem_u_new, std::vector<cipher*>& mem_y_pre, std::vector<cipher*>& mem_u_pre)
        {
            for(int i = 0; i < 8; i++)
            {
                if(i < 4)
                {
                    this->mem_stack_new[2 * i]->coeff = mem_y_new[i]->ciphertext[0]->coeff;
                    this->mem_stack_new[2 * i + 1]->coeff = mem_y_new[i]->ciphertext[1]->coeff;

                    this->mem_stack_pre[2 * i]->coeff = mem_y_pre[i]->ciphertext[0]->coeff;
                    this->mem_stack_pre[2 * i + 1]->coeff = mem_y_pre[i]->ciphertext[1]->coeff;
                }
                else
                {
                    this->mem_stack_new[2 * i]->coeff = mem_u_new[i - 4]->ciphertext[0]->coeff;
                    this->mem_stack_new[2 * i + 1]->coeff = mem_u_new[i - 4]->ciphertext[1]->coeff;

                    this->mem_stack_pre[2 * i]->coeff = mem_u_pre[i - 4]->ciphertext[0]->coeff;
                    this->mem_stack_pre[2 * i + 1]->coeff = mem_u_pre[i - 4]->ciphertext[1]->coeff;
                }
            }

            if(this->timing)
            { 
                // this is odd section
                poly_handler::poly_cascade_concat(this->mem_stack_new, this->arx_state_even);
                poly_handler::poly_cascade_concat(this->mem_stack_pre, this->arx_state_odd);
                group_handler::group_dot(this->g_r_1, arx_state_even, this->odd_nu[0]);
                group_handler::group_dot(this->g_gamma_1_r_1, arx_state_even, this->odd_nu[1]);
                group_handler::group_dot(this->g_alpha_1_r_1_F_r_0_m, arx_state_odd, this->odd_nu[2]);
                group_handler::group_dot(this->g_r_1_F_r_0_m, arx_state_odd, this->odd_nu[3]);

                group_handler::group_dot(this->g_beta_1_s_H_r_0_m, arx_state_odd, this->odd_mu[0]);
                group_handler::group_dot(this->g_s_H_r_0_m, arx_state_odd, this->odd_mu[1]);
            }
            else
            {
                // this is even section
                poly_handler::poly_cascade_concat(this->mem_stack_new, this->arx_state_odd);
                poly_handler::poly_cascade_concat(this->mem_stack_pre, this->arx_state_even);
                group_handler::group_dot(this->g_r_0, arx_state_odd, this->even_nu[0]);
                group_handler::group_dot(this->g_gamma_0_r_0, arx_state_odd, this->even_nu[1]);
                group_handler::group_dot(this->g_alpha_0_r_0_F_r_1_m, arx_state_even, this->even_nu[2]);
                group_handler::group_dot(this->g_r_0_F_r_1_m, arx_state_even, this->even_nu[3]);
                group_handler::group_dot(this->g_beta_0_s_H_r_1_m, arx_state_even, this->even_mu[0]);
                group_handler::group_dot(this->g_s_H_r_1_m, arx_state_even, this->even_mu[1]);
            }
        }

        bool verifying_proof(cipher* y, cipher* u_re_enc, cipher* u, mpz_class& previous_pf)
        {
            poly_handler::poly_cascade_concat(y->ciphertext, this->Y);
            poly_handler::poly_cascade_concat(u_re_enc->ciphertext, this->U_re);
            poly_handler::poly_cascade_concat(u->ciphertext, this->U);
            
            bool check = false;

            if(this->timing)
            {
                // this is odd section
                mpz_class vc1, vc2, vc3;
                mpz_class pf1, pf2, pf3, pf4, pf5;
                
                group_handler::group_dot(this->g_r_1_G, Y, vc1);
                group_handler::group_dot(this->g_r_1_R, U_re, vc2);
                group_handler::group_dot(this->g_s, U, vc3);
                vc3 = ((vc3) % this->g_mod + this->g_mod) % this->g_mod;

                mpz_powm(pf1.get_mpz_t(), this->odd_nu[0].get_mpz_t(), this->gamma_1->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                mpz_powm(pf2.get_mpz_t(), this->odd_nu[3].get_mpz_t(), this->alpha_1->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                mpz_powm(pf3.get_mpz_t(), this->odd_mu[1].get_mpz_t(), this->beta_1->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                pf4 = (this->odd_nu[3] * previous_pf) % this->g_mod;
                pf4 = (pf4 * vc1) % this->g_mod;
                pf4 = ((pf4 * vc2) % this->g_mod + this->g_mod) % this->g_mod;
                pf5 = ((this->odd_mu[1] * previous_pf) % this->g_mod + this->g_mod) % this->g_mod;

                if(
                    pf1 != this->odd_nu[1] ||
                    pf2 != this->odd_nu[2] ||
                    pf3 != this->odd_mu[0] ||
                    pf4 != this->odd_nu[0] ||
                    pf5 != vc3
                )
                {
                    std::cout << " ------------------------------- " << std::endl;
                    std::cout << "  x(t+1) linearity check : " << pf1 - this->odd_nu[1] << std::endl;
                    std::cout << "  (rF - r) linearity check : " << pf2 - this->odd_nu[2] << std::endl;
                    std::cout << "  (sH - r) linearity check : " << pf3 - this->odd_mu[0] << std::endl;
                    std::cout << "  x(t+1) = Fx(t) + Gy + Ru check : " << pf4 - this->odd_nu[0] << std::endl;
                    std::cout << "  u(t) = Hx(t) check : " << pf5 - vc3 << std::endl;
                    std::cout << " ------------------------------- " << std::endl;
                    check = false;
                }
                else
                {
                    check = true;
                }
            }
            else
            {
                // this is even section
                mpz_class vc1, vc2, vc3;
                mpz_class pf1, pf2, pf3, pf4, pf5;

                group_handler::group_dot(this->g_r_0_G, Y, vc1);
                group_handler::group_dot(this->g_r_0_R, U_re, vc2);
                group_handler::group_dot(this->g_s, U, vc3);
                vc3 = ((vc3) % this->g_mod + this->g_mod) % this->g_mod;

                mpz_powm(pf1.get_mpz_t(), this->even_nu[0].get_mpz_t(), this->gamma_0->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                mpz_powm(pf2.get_mpz_t(), this->even_nu[3].get_mpz_t(), this->alpha_0->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                mpz_powm(pf3.get_mpz_t(), this->even_mu[1].get_mpz_t(), this->beta_0->coeff[0].get_mpz_t(), this->g_mod.get_mpz_t());
                pf4 = (this->even_nu[3] * previous_pf) % this->g_mod;
                pf4 = (pf4 * vc1) % this->g_mod;
                pf4 = ((pf4 * vc2) % this->g_mod + this->g_mod) % this->g_mod;
                pf5 = ((this->even_mu[1] * previous_pf) % this->g_mod + this->g_mod) % this->g_mod;

                if(
                    pf1 != this->even_nu[1] ||
                    pf2 != this->even_nu[2] ||
                    pf3 != this->even_mu[0] ||
                    pf4 != this->even_nu[0] ||
                    pf5 != vc3
                )
                {
                    std::cout << " ------------------------------- " << std::endl;
                    std::cout << "  x(t+1) linearity check : " << pf1 - this->even_nu[1] << std::endl;
                    std::cout << "  (rF - r) linearity check : " << pf2 - this->even_nu[2] << std::endl;
                    std::cout << "  (sH - r) linearity check : " << pf3 - this->even_mu[0] << std::endl;
                    std::cout << "  x(t+1) = Fx(t) + Gy + Ru check : " << pf4 - this->even_nu[0] << std::endl;
                    std::cout << "  u(t) = Hx(t) check : " << pf5 - vc3 << std::endl;
                    std::cout << " ------------------------------- " << std::endl;
                    check = false;
                }
                else
                {
                    check = true;
                }
            }

            previous_pf = this->timing ? this->odd_nu[0] : this->even_nu[0];
            this->timing = !this->timing;
            return check;
        }
};

class authentic_ecc
{
    public:
        int poly_degree;
        mpz_class cipher_mod;
        mpz_class ecc_mod;

        // ekf contain(on ECC) gamma_0*r_0, alpha_0(r_0*F - r_1), (r_0*F - r_1), r_0*G, r_0*R,
        //                       gamma_1*r_1, alpha_1(r_1*F - r_0), (r_1*F - r_0), r_1*G, r_1*R,
        //                       beta_0(s*H - r_1), (s*H - r_1), beta_1(s*H - r_0), (s*H - r_0),
        //                       r_0, r_1, s, respectively. (only memory management)
        std::vector<eccvec*> ekf;
        // just pointer mapping (to use easy)
        eccvec* ecc_gamma_0_r_0;
        eccvec* ecc_alpha_0_r_0_F_r_1_m;
        eccvec* ecc_r_0_F_r_1_m;
        eccvec* ecc_r_0_G;
        eccvec* ecc_r_0_R;
        eccvec* ecc_gamma_1_r_1;
        eccvec* ecc_alpha_1_r_1_F_r_0_m;
        eccvec* ecc_r_1_F_r_0_m;
        eccvec* ecc_r_1_G;
        eccvec* ecc_r_1_R;
        eccvec* ecc_beta_0_s_H_r_1_m;
        eccvec* ecc_s_H_r_1_m;
        eccvec* ecc_beta_1_s_H_r_0_m;
        eccvec* ecc_s_H_r_0_m;
        eccvec* ecc_r_0;
        eccvec* ecc_r_1;
        eccvec* ecc_s;

        // key contain r_0, r_1, s, Gamma_0, Gamma_1, alpha_0, alpha_1, beta_0, beta_1, respectively. (only memory management)
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


        std::vector<poly*> mem_stack_new;
        std::vector<poly*> mem_stack_pre;

        poly* arx_state_even;
        poly* arx_state_odd;

        bool timing = false;
        std::vector<point> even_nu;
        std::vector<point> even_mu;
        std::vector<point> odd_nu;
        std::vector<point> odd_mu;

        poly* Y;
        poly* U_re;
        poly* U;

        authentic_ecc(int poly_degree, mpz_class cipher_mod):
            poly_degree(poly_degree),
            cipher_mod(cipher_mod)
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

            this->mem_stack_new.resize(16);
            this->mem_stack_pre.resize(16);

            poly* temp;
            for(int i = 0; i < 16; i++)
            {
                temp = new poly(this->poly_degree);
                this->mem_stack_new[i] = temp;
                temp = new poly(this->poly_degree);
                this->mem_stack_pre[i] = temp;
            }

            this->arx_state_even = new poly(this->poly_degree * this->mem_stack_new.size());
            this->arx_state_odd = new poly(this->poly_degree * this->mem_stack_new.size());

            this->even_nu.resize(4);
            this->even_mu.resize(2);
            this->odd_nu.resize(4);
            this->odd_mu.resize(2);

            this->Y = new poly(this->poly_degree * 2);
            this->U_re = new poly(this->poly_degree * 2);
            this->U = new poly(this->poly_degree * 3);
        };

        ~authentic_ecc()
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

            for(auto& factor : this->mem_stack_new)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_stack_new.clear();
            for(auto& factor : this->mem_stack_pre)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_stack_pre.clear();

            delete this->arx_state_even;
            delete this->arx_state_odd;

            this->even_nu.clear();
            this->even_mu.clear();
            this->odd_nu.clear();
            this->odd_mu.clear();

            delete this->Y;
            delete this->U_re;
            delete this->U;
        }

        void make_ekf(std::vector<cipher*>& P_enc, std::vector<cipher*>& Q_enc)
        {
            // ciphertext whole size(N size message part and N size random key part, respectivley)
            int step_size = 2 * this->poly_degree;

            // for transform to eccvec
            eccvec* table = new eccvec(256); //bit size
            ecc_handler::table_set(table);

            // for save ecc mod
            this->ecc_mod = table->ecc_mod;

            // gamma0 * r0 power g
            poly* temp_gamma_0_r_0 = new poly(this->r_0->ring_dim);
            poly_handler::poly_scal_mul_p(this->gamma_0, this->r_0, temp_gamma_0_r_0);
            poly_handler::poly_mod(temp_gamma_0_r_0, this->cipher_mod, temp_gamma_0_r_0);
            eccvec* gamma_0_r_0 = new eccvec(this->r_0->ring_dim);
            ecc_handler::poly_2_ecc(temp_gamma_0_r_0, table, gamma_0_r_0);
            this->ecc_gamma_0_r_0 = gamma_0_r_0;
            this->ekf.push_back(gamma_0_r_0);
            delete temp_gamma_0_r_0;

            // alpha0 * (r0 * F - r1) power g and (r0 * F - r1) power g
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
            poly_handler::poly_mod(temp_alpha_0_r_0_F_r_1_m, this->cipher_mod, temp_alpha_0_r_0_F_r_1_m);
            eccvec* r_0_F_r_1_m = new eccvec(this->r_0->ring_dim);
            ecc_handler::poly_2_ecc(temp_alpha_0_r_0_F_r_1_m, table, r_0_F_r_1_m);
            this->ecc_r_0_F_r_1_m = r_0_F_r_1_m;
            poly_handler::poly_scal_mul_p(this->alpha_0, temp_alpha_0_r_0_F_r_1_m, temp_alpha_0_r_0_F_r_1_m);
            poly_handler::poly_mod(temp_alpha_0_r_0_F_r_1_m, this->cipher_mod, temp_alpha_0_r_0_F_r_1_m);
            eccvec* alpha_0_r_0_F_r_1_m = new eccvec(this->r_0->ring_dim);
            ecc_handler::poly_2_ecc(temp_alpha_0_r_0_F_r_1_m, table, alpha_0_r_0_F_r_1_m);
            this->ecc_alpha_0_r_0_F_r_1_m = alpha_0_r_0_F_r_1_m;
            this->ekf.push_back(alpha_0_r_0_F_r_1_m);
            this->ekf.push_back(r_0_F_r_1_m);
            delete temp_alpha_0_r_0_F_r_1_m;
            delete temp_r_1_m1;

            // r_0 * G power g
            poly* temp_r_0_G = new poly(step_size);
            temp_r_0_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_0_G->coeff[i] = this->r_0->coeff[i + 3 * step_size];
            }
            eccvec* r_0_G = new eccvec(step_size);
            ecc_handler::poly_2_ecc(temp_r_0_G, table, r_0_G);
            this->ecc_r_0_G = r_0_G;
            this->ekf.push_back(r_0_G);
            delete temp_r_0_G;

            // r_0 * R power g
            poly* temp_r_0_R = new poly(step_size);
            temp_r_0_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_0_R->coeff[i] = this->r_0->coeff[i + 7 * step_size];
            }
            eccvec* r_0_R = new eccvec(step_size);
            ecc_handler::poly_2_ecc(temp_r_0_R, table, r_0_R);
            this->ecc_r_0_R = r_0_R;
            this->ekf.push_back(r_0_R);
            delete temp_r_0_R;

            // gamma_1 * r_1 power g
            poly* temp_gamma_1_r_1 = new poly(this->r_1->ring_dim);
            poly_handler::poly_scal_mul_p(this->gamma_1, this->r_1, temp_gamma_1_r_1);
            poly_handler::poly_mod(temp_gamma_1_r_1, this->cipher_mod, temp_gamma_1_r_1);
            eccvec* gamma_1_r_1 = new eccvec(this->r_1->ring_dim);
            ecc_handler::poly_2_ecc(temp_gamma_1_r_1, table, gamma_1_r_1);
            this->ecc_gamma_1_r_1 = gamma_1_r_1;
            this->ekf.push_back(gamma_1_r_1);
            delete temp_gamma_1_r_1;

            // alpha_1 * (r_1 * F - r_0) power g and (r_0 * F - r_1) power g
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
            poly_handler::poly_mod(temp_alpha_1_r_1_F_r_0_m, this->cipher_mod, temp_alpha_1_r_1_F_r_0_m);
            eccvec* r_1_F_r_0_m = new eccvec(this->r_1->ring_dim);
            ecc_handler::poly_2_ecc(temp_alpha_1_r_1_F_r_0_m, table, r_1_F_r_0_m);
            this->ecc_r_1_F_r_0_m = r_1_F_r_0_m;
            poly_handler::poly_scal_mul_p(this->alpha_1, temp_alpha_1_r_1_F_r_0_m, temp_alpha_1_r_1_F_r_0_m);
            poly_handler::poly_mod(temp_alpha_1_r_1_F_r_0_m, this->cipher_mod, temp_alpha_1_r_1_F_r_0_m);
            eccvec* alpha_1_r_1_F_r_0_m = new eccvec(this->r_1->ring_dim);
            ecc_handler::poly_2_ecc(temp_alpha_1_r_1_F_r_0_m, table, alpha_1_r_1_F_r_0_m);
            this->ecc_alpha_1_r_1_F_r_0_m = alpha_1_r_1_F_r_0_m;
            this->ekf.push_back(alpha_1_r_1_F_r_0_m);
            this->ekf.push_back(r_1_F_r_0_m);
            delete temp_alpha_1_r_1_F_r_0_m;
            delete temp_r_0_m1;

            // r_1 * G power g
            poly* temp_r_1_G = new poly(step_size);
            temp_r_1_G->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_1_G->coeff[i] = this->r_1->coeff[i + 3 * step_size];
            }
            eccvec* r_1_G = new eccvec(step_size);
            ecc_handler::poly_2_ecc(temp_r_1_G, table, r_1_G);
            this->ecc_r_1_G = r_1_G;
            this->ekf.push_back(r_1_G);
            delete temp_r_1_G;

            // r_1 * R power g
            poly* temp_r_1_R = new poly(step_size);
            temp_r_1_R->fill_zero();
            for(int i = 0; i < step_size; i++)
            {
                temp_r_1_R->coeff[i] = this->r_1->coeff[i + 7 * step_size];
            }
            eccvec* r_1_R = new eccvec(step_size);
            ecc_handler::poly_2_ecc(temp_r_1_R, table, r_1_R);
            this->ecc_r_1_R = r_1_R;
            this->ekf.push_back(r_1_R);
            delete temp_r_1_R;

            // beta_0(s*H - r_1), (s*H - r_1), beta_1(s*H - r_0), (s*H - r_0), respectively
            // mr_cipher* comp_H_enc_mat =  new mr_cipher(P_enc[0]->ring_dim, P_enc[0]->plain_mod, P_enc[0]->cipher_mod);
            std::vector<poly*> temp_sH(8);
            for(int i = 0; i < 8; i++)
            {
                temp_sH[i] = new poly(2 * this->poly_degree);
                if(i < 4)
                {
                    // format_transform_handler::cipher_2_mr_cipher(P_enc[i], comp_H_enc_mat);
                    // crypto_handler::pval_mr_mul(this->s, comp_H_enc_mat, temp_sH[i]);
                    crypto_handler::pval_mrlike_mul(this->s, P_enc[i], temp_sH[i]);
                }
                else
                {
                    // format_transform_handler::cipher_2_mr_cipher(Q_enc[i - 4], comp_H_enc_mat);
                    // crypto_handler::pval_mr_mul(this->s, comp_H_enc_mat, temp_sH[i]);
                    crypto_handler::pval_mrlike_mul(this->s, Q_enc[i - 4], temp_sH[i]);
                }
            }

            // delete comp_H_enc_mat;

            poly* sH = poly_handler::poly_recur_concat(temp_sH);

            for(auto& factor : temp_sH)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            temp_sH.clear();

            poly* temp_r_1_m2 = new poly(this->r_1->ring_dim);
            poly_handler::poly_neg(this->r_1, temp_r_1_m2);
            poly* temp_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_1_m2, temp_s_H_r_1_m);
            poly_handler::poly_mod(temp_s_H_r_1_m, this->cipher_mod, temp_s_H_r_1_m);
            eccvec* s_H_r_1_m = new eccvec(sH->ring_dim);
            ecc_handler::poly_2_ecc(temp_s_H_r_1_m, table, s_H_r_1_m);
            this->ecc_s_H_r_1_m = s_H_r_1_m;
            poly* temp_beta_0_s_H_r_1_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_0, temp_s_H_r_1_m, temp_beta_0_s_H_r_1_m);
            poly_handler::poly_mod(temp_beta_0_s_H_r_1_m, this->cipher_mod, temp_beta_0_s_H_r_1_m);
            eccvec* beta_0_s_H_r_1_m = new eccvec(sH->ring_dim);
            ecc_handler::poly_2_ecc(temp_beta_0_s_H_r_1_m, table, beta_0_s_H_r_1_m);
            this->ecc_beta_0_s_H_r_1_m = beta_0_s_H_r_1_m;
            this->ekf.push_back(beta_0_s_H_r_1_m);
            this->ekf.push_back(s_H_r_1_m);
            delete temp_r_1_m2;
            delete temp_s_H_r_1_m;
            delete temp_beta_0_s_H_r_1_m;
            
            poly* temp_r_0_m2 = new poly(this->r_0->ring_dim);
            poly_handler::poly_neg(this->r_0, temp_r_0_m2);
            poly* temp_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_add(sH, temp_r_0_m2, temp_s_H_r_0_m);
            poly_handler::poly_mod(temp_s_H_r_0_m, this->cipher_mod, temp_s_H_r_0_m);
            eccvec* s_H_r_0_m = new eccvec(sH->ring_dim);
            ecc_handler::poly_2_ecc(temp_s_H_r_0_m, table, s_H_r_0_m);
            this->ecc_s_H_r_0_m = s_H_r_0_m;
            poly* temp_beta_1_s_H_r_0_m = new poly(sH->ring_dim);
            poly_handler::poly_scal_mul_p(this->beta_1, temp_s_H_r_0_m, temp_beta_1_s_H_r_0_m);
            poly_handler::poly_mod(temp_beta_1_s_H_r_0_m, this->cipher_mod, temp_beta_1_s_H_r_0_m);
            eccvec* beta_1_s_H_r_0_m = new eccvec(sH->ring_dim);
            ecc_handler::poly_2_ecc(temp_beta_1_s_H_r_0_m, table, beta_1_s_H_r_0_m);
            this->ecc_beta_1_s_H_r_0_m = beta_1_s_H_r_0_m;
            this->ekf.push_back(beta_1_s_H_r_0_m);
            this->ekf.push_back(s_H_r_0_m);
            delete temp_r_0_m2;
            delete temp_s_H_r_0_m;
            delete temp_beta_1_s_H_r_0_m;
            delete sH;
            
            // r_0 power g
            eccvec* r_0 = new eccvec(this->r_0->ring_dim);
            ecc_handler::poly_2_ecc(this->r_0, table, r_0);
            this->ecc_r_0 = r_0;
            this->ekf.push_back(r_0);

            // r_1 power g
            eccvec* r_1 = new eccvec(this->r_1->ring_dim);
            ecc_handler::poly_2_ecc(this->r_1, table, r_1);
            this->ecc_r_1 = r_1;
            this->ekf.push_back(r_1);

            // s power g
            eccvec* s = new eccvec(this->s->ring_dim);
            ecc_handler::poly_2_ecc(this->s, table, s);
            this->ecc_s = s;
            this->ekf.push_back(s);

            delete table;
        }

        void generate_proof(std::vector<cipher*>& mem_y_new, std::vector<cipher*>& mem_u_new, std::vector<cipher*>& mem_y_pre, std::vector<cipher*>& mem_u_pre)
        {
            for(int i = 0; i < 8; i++)
            {
                if(i < 4)
                {
                    this->mem_stack_new[2 * i]->coeff = mem_y_new[i]->ciphertext[0]->coeff;
                    this->mem_stack_new[2 * i + 1]->coeff = mem_y_new[i]->ciphertext[1]->coeff;

                    this->mem_stack_pre[2 * i]->coeff = mem_y_pre[i]->ciphertext[0]->coeff;
                    this->mem_stack_pre[2 * i + 1]->coeff = mem_y_pre[i]->ciphertext[1]->coeff;
                }
                else
                {
                    this->mem_stack_new[2 * i]->coeff = mem_u_new[i - 4]->ciphertext[0]->coeff;
                    this->mem_stack_new[2 * i + 1]->coeff = mem_u_new[i - 4]->ciphertext[1]->coeff;

                    this->mem_stack_pre[2 * i]->coeff = mem_u_pre[i - 4]->ciphertext[0]->coeff;
                    this->mem_stack_pre[2 * i + 1]->coeff = mem_u_pre[i - 4]->ciphertext[1]->coeff;
                }
            }

            if(this->timing)
            { 
                // this is odd section
                poly_handler::poly_cascade_concat(this->mem_stack_new, this->arx_state_even);
                poly_handler::poly_cascade_concat(this->mem_stack_pre, this->arx_state_odd);

                // ecc_handler::ecc_dot(this->ecc_r_1, arx_state_even, this->odd_nu[0]);
                // ecc_handler::ecc_dot(this->ecc_gamma_1_r_1, arx_state_even, this->odd_nu[1]);
                // ecc_handler::ecc_dot(this->ecc_alpha_1_r_1_F_r_0_m, arx_state_odd, this->odd_nu[2]);
                // ecc_handler::ecc_dot(this->ecc_r_1_F_r_0_m, arx_state_odd, this->odd_nu[3]);
                // ecc_handler::ecc_dot(this->ecc_beta_1_s_H_r_0_m, arx_state_odd, this->odd_mu[0]);
                // ecc_handler::ecc_dot(this->ecc_s_H_r_0_m, arx_state_odd, this->odd_mu[1]);

                ecc_handler::ecc_dot_j(this->ecc_r_1, arx_state_even, this->odd_nu[0]);
                ecc_handler::ecc_dot_j(this->ecc_gamma_1_r_1, arx_state_even, this->odd_nu[1]);
                ecc_handler::ecc_dot_j(this->ecc_alpha_1_r_1_F_r_0_m, arx_state_odd, this->odd_nu[2]);
                ecc_handler::ecc_dot_j(this->ecc_r_1_F_r_0_m, arx_state_odd, this->odd_nu[3]);
                ecc_handler::ecc_dot_j(this->ecc_beta_1_s_H_r_0_m, arx_state_odd, this->odd_mu[0]);
                ecc_handler::ecc_dot_j(this->ecc_s_H_r_0_m, arx_state_odd, this->odd_mu[1]);
            }
            else
            {
                // this is even section
                poly_handler::poly_cascade_concat(this->mem_stack_new, this->arx_state_odd);
                poly_handler::poly_cascade_concat(this->mem_stack_pre, this->arx_state_even);

                // ecc_handler::ecc_dot(this->ecc_r_0, arx_state_odd, this->even_nu[0]);
                // ecc_handler::ecc_dot(this->ecc_gamma_0_r_0, arx_state_odd, this->even_nu[1]);
                // ecc_handler::ecc_dot(this->ecc_alpha_0_r_0_F_r_1_m, arx_state_even, this->even_nu[2]);
                // ecc_handler::ecc_dot(this->ecc_r_0_F_r_1_m, arx_state_even, this->even_nu[3]);
                // ecc_handler::ecc_dot(this->ecc_beta_0_s_H_r_1_m, arx_state_even, this->even_mu[0]);
                // ecc_handler::ecc_dot(this->ecc_s_H_r_1_m, arx_state_even, this->even_mu[1]);

                ecc_handler::ecc_dot_j(this->ecc_r_0, arx_state_odd, this->even_nu[0]);
                ecc_handler::ecc_dot_j(this->ecc_gamma_0_r_0, arx_state_odd, this->even_nu[1]);
                ecc_handler::ecc_dot_j(this->ecc_alpha_0_r_0_F_r_1_m, arx_state_even, this->even_nu[2]);
                ecc_handler::ecc_dot_j(this->ecc_r_0_F_r_1_m, arx_state_even, this->even_nu[3]);
                ecc_handler::ecc_dot_j(this->ecc_beta_0_s_H_r_1_m, arx_state_even, this->even_mu[0]);
                ecc_handler::ecc_dot_j(this->ecc_s_H_r_1_m, arx_state_even, this->even_mu[1]);
            }
        }

        bool verifying_proof(cipher* y, cipher* u_re_enc, cipher* u, point& previous_pf)
        {
            poly_handler::poly_cascade_concat(y->ciphertext, this->Y);
            poly_handler::poly_cascade_concat(u_re_enc->ciphertext, this->U_re);
            poly_handler::poly_cascade_concat(u->ciphertext, this->U);
            
            bool check = false;

            if(this->timing)
            {
                // this is odd section
                point vc1, vc2, vc3;
                point pf1, pf2, pf3, pf4, pf5;
                
                // ecc_handler::ecc_dot(this->ecc_r_1_G, Y, vc1);
                // ecc_handler::ecc_dot(this->ecc_r_1_R, U_re, vc2);
                // ecc_handler::ecc_dot(this->ecc_s, U, vc3); 

                ecc_handler::ecc_dot_j(this->ecc_r_1_G, Y, vc1);
                ecc_handler::ecc_dot_j(this->ecc_r_1_R, U_re, vc2);
                ecc_handler::ecc_dot_j(this->ecc_s, U, vc3); 

                ecc_handler::point_mul(this->odd_nu[0], this->gamma_1->coeff[0], pf1, this->ecc_mod);
                ecc_handler::point_mul(this->odd_nu[3], this->alpha_1->coeff[0], pf2, this->ecc_mod);
                ecc_handler::point_mul(this->odd_mu[1], this->beta_1->coeff[0], pf3, this->ecc_mod);
                ecc_handler::point_add(this->odd_nu[3], previous_pf, pf4, this->ecc_mod);
                ecc_handler::point_add(pf4, vc1, pf4, this->ecc_mod);
                ecc_handler::point_add(pf4, vc2, pf4, this->ecc_mod);
                ecc_handler::point_add(this->odd_mu[1], previous_pf, pf5, this->ecc_mod);

                if(
                    !ecc_handler::compare(pf1, this->odd_nu[1]) ||
                    !ecc_handler::compare(pf2, this->odd_nu[2]) ||
                    !ecc_handler::compare(pf3, this->odd_mu[0]) ||
                    !ecc_handler::compare(pf4, this->odd_nu[0]) ||
                    !ecc_handler::compare(pf5, vc3) 
                )
                {
                    std::cout << " ------------------------------- " << std::endl;
                    std::cout << "  x(t+1) linearity check : " << !ecc_handler::compare(pf1, this->odd_nu[1]) << std::endl;
                    std::cout << "  (rF - r) linearity check : " << !ecc_handler::compare(pf2, this->odd_nu[2]) << std::endl;
                    std::cout << "  (sH - r) linearity check : " << !ecc_handler::compare(pf3, this->odd_mu[0]) << std::endl;
                    std::cout << "  x(t+1) = Fx(t) + Gy + Ru check : " << !ecc_handler::compare(pf4, this->odd_nu[0]) << std::endl;
                    std::cout << "  u(t) = Hx(t) check : " << !ecc_handler::compare(pf5, vc3) << std::endl;
                    std::cout << " ------------------------------- " << std::endl;
                    check = false;
                }
                else
                {
                    check = true;
                }
            }
            else
            {
                // this is even section
                point vc1, vc2, vc3;
                point pf1, pf2, pf3, pf4, pf5;

                // ecc_handler::ecc_dot(this->ecc_r_0_G, Y, vc1);
                // ecc_handler::ecc_dot(this->ecc_r_0_R, U_re, vc2);
                // ecc_handler::ecc_dot(this->ecc_s, U, vc3);

                ecc_handler::ecc_dot_j(this->ecc_r_0_G, Y, vc1);
                ecc_handler::ecc_dot_j(this->ecc_r_0_R, U_re, vc2);
                ecc_handler::ecc_dot_j(this->ecc_s, U, vc3);

                ecc_handler::point_mul(this->even_nu[0], this->gamma_0->coeff[0], pf1, this->ecc_mod);
                ecc_handler::point_mul(this->even_nu[3], this->alpha_0->coeff[0], pf2, this->ecc_mod);
                ecc_handler::point_mul(this->even_mu[1], this->beta_0->coeff[0], pf3, this->ecc_mod);
                ecc_handler::point_add(this->even_nu[3], previous_pf, pf4, this->ecc_mod);
                ecc_handler::point_add(pf4, vc1, pf4, this->ecc_mod);
                ecc_handler::point_add(pf4, vc2, pf4, this->ecc_mod);
                ecc_handler::point_add(this->even_mu[1], previous_pf, pf5, this->ecc_mod);

                if(
                    !ecc_handler::compare(pf1, this->even_nu[1]) ||
                    !ecc_handler::compare(pf2, this->even_nu[2]) ||
                    !ecc_handler::compare(pf3, this->even_mu[0]) ||
                    !ecc_handler::compare(pf4, this->even_nu[0]) ||
                    !ecc_handler::compare(pf5, vc3) 
                )
                {
                    std::cout << " ------------------------------- " << std::endl;
                    std::cout << "  x(t+1) linearity check : " << !ecc_handler::compare(pf1, this->even_nu[1]) << std::endl;
                    std::cout << "  (rF - r) linearity check : " << !ecc_handler::compare(pf2, this->even_nu[2]) << std::endl;
                    std::cout << "  (sH - r) linearity check : " << !ecc_handler::compare(pf3, this->even_mu[0]) << std::endl;
                    std::cout << "  x(t+1) = Fx(t) + Gy + Ru check : " << !ecc_handler::compare(pf4, this->even_nu[0]) << std::endl;
                    std::cout << "  u(t) = Hx(t) check : " << !ecc_handler::compare(pf5, vc3) << std::endl;
                    std::cout << " ------------------------------- " << std::endl;
                    check = false;
                }
                else
                {
                    check = true;
                }
            }

            previous_pf = this->timing ? this->odd_nu[0] : this->even_nu[0];
            this->timing = !this->timing;
            return check;
        }
};

#endif