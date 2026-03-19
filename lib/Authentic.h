#ifndef AUTHENTIC_H
#define AUTHENTIC_H

#include "Group.h"
#include "RLWE.h"

class authentic
{
    private: 
        mpz_class cipher_mod;
        mpz_class g_mod;
        mpz_class g_gen;
        
        std::vector<gvec*> ekf;
        std::vector<poly*> key;
        
    public:
        authentic(int poly_degree, mpz_class cipher_mod, mpz_class g_mod, mpz_class g_gen):
            g_mod(g_mod),
            g_gen(g_gen)
        {
            poly* r_0 = new poly(24 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, r_0);
            this->key.push_back(r_0);

            poly* r_1 = new poly(24 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, r_1);
            this->key.push_back(r_1);

            poly* s = new poly(3 * poly_degree);
            random_handler::not_zero_public_key(cipher_mod, s);
            this->key.push_back(s);

            poly* alpha_0 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, alpha_0);
            this->key.push_back(alpha_0);

            poly* alpha_1 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, alpha_1);
            this->key.push_back(alpha_1);

            poly* beta_0 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, beta_0);
            this->key.push_back(beta_0);

            poly* beta_1 = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, beta_1);
            this->key.push_back(beta_1);

            poly* delta = new poly(1);
            random_handler::not_zero_public_key(cipher_mod, delta);
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

        void make_ekf()
        {}
};

#endif
