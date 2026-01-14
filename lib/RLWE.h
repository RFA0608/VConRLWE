#ifndef RLWE_H
#define RLWE_H

#include <random>
#include <ctime>

#include "Struct.h"

class cipher
{
    public:
        std::vector<poly*> ciphertext;
        int ring_dim;
        int level = 0;
        mpz_class plain_mod;
        mpz_class cipher_mod;
        mpz_class psi_plain;
        mpz_class psi_cipher;
        // mpz_class g_plain;
        // mpz_class g_cipher;

        cipher(int ring_dim, mpz_class plain_mod, mpz_class cipher_mod, mpz_class psi_p, mpz_class psi_c):
            ring_dim(ring_dim), 
            plain_mod(plain_mod),
            cipher_mod(cipher_mod),
            psi_plain(psi_p),
            psi_cipher(psi_c)
        {}

        ~cipher()
        {
            for(auto& factor : this->ciphertext)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->ciphertext.clear();
        }

        cipher* clone()
        {
            cipher* clone = new cipher(this->ring_dim, this->plain_mod, this->cipher_mod, this->psi_plain, this->psi_cipher);

            clone->level = this->level;
            for(auto& factor : this->ciphertext)
            {
                if(factor != nullptr)
                {
                    clone->ciphertext.push_back(factor->clone());
                }
            }

            return clone;
        }
};

class random_handler
{
    public:
        static void sample_cbd(std::mt19937& gen, mpz_class& res)
        {
            std::uniform_int_distribution<int> dist(0, 1);

            int a = dist(gen) + dist(gen);
            int b = dist(gen) + dist(gen);

            res = a - b;
        }

        static void sample_uni(const mpz_class& mod, gmp_randstate_t& state, mpz_class& res)
        {
            mpz_urandomm(res.get_mpz_t(), state, mod.get_mpz_t());
        }

        static void secret_key(poly* sk)
        {
            thread_local std::mt19937 gen(std::random_device{}());

            mpz_class temp;

            for(int i = 0; i < sk->ring_dim; i++)
            {
                sample_cbd(gen, temp);
                sk->coeff[i] = temp;
            }
        }

        static void public_key(const mpz_class& mod, poly* pb)
        {
            mpz_class temp;

            for(int i = 0; i < pb->ring_dim; i++)
            {
                sample_uni(mod, seed_gen::get_rand_state(), temp);
                pb->coeff[i] = temp;
            }
        }

        static void errors(poly* err)
        {
            thread_local std::mt19937 gen(std::random_device{}());

            mpz_class temp;

            for(int i = 0; i < err->ring_dim; i++)
            {
                sample_cbd(gen, temp);
                err->coeff[i] = temp;
            }
        }
};

class crypto_handler
{
    public:
        static int encrypt(poly* m, poly* sk, cipher* res)
        {
            if((res->ring_dim != m->ring_dim) || (res->ring_dim != sk->ring_dim))
            {
                return -1;
            }
            else
            {

                switch(res->level)
                {
                    case 0:
                    {
                        res->level = 1;

                        break;
                    }
                    case 1:
                    {
                        delete res->ciphertext.back();
                        res->ciphertext.pop_back();
                        delete res->ciphertext.back();
                        res->ciphertext.pop_back();

                        break;
                    }
                    default:
                    {
                        return -1;
                    }
                }

                poly* a = new poly(res->ring_dim);
                random_handler::public_key(res->cipher_mod, a);

                poly* ska = new poly(res->ring_dim);
                poly_handler::poly_mul_ntt(sk, a, res->cipher_mod, res->psi_cipher, ska);

                poly* e = new poly(res->ring_dim);
                random_handler::errors(e);
                poly_handler::poly_scal_mul(res->plain_mod, e, e);

                poly* b = new poly(res->ring_dim);
                poly_handler::poly_add(m, ska, b);
                poly_handler::poly_add(b, e, b);
                poly_handler::poly_mod(b, res->cipher_mod, b);

                res->ciphertext.push_back(b);
                res->ciphertext.push_back(a);
                
                delete ska;
                delete e;

                return 0;
            }
        }

        static int decrypt(cipher* c, poly* sk, poly* res)
        {
            if((res->ring_dim != c->ring_dim) || (res->ring_dim != sk->ring_dim))
            {
                return -1;
            }
            else
            {
                switch(c->level)
                {
                    case 1:
                    {
                        poly* ska = new poly(res->ring_dim);

                        poly_handler::poly_mul_ntt(sk, c->ciphertext[1], c->cipher_mod, c->psi_cipher, ska);
                        poly_handler::poly_neg(ska, ska);
                        poly_handler::poly_add(c->ciphertext[0], ska, res);
                        poly_handler::poly_mod(res, c->cipher_mod, res);
                        poly_handler::poly_bias_mod(res, c->cipher_mod, res);
                        poly_handler::poly_mod(res, c->plain_mod, res);
                        
                        delete ska;
                        break;
                    }
                    case 2:
                    {
                        poly* sk2 = new poly(res->ring_dim);
                        poly_handler::poly_mul_ntt(sk, sk, c->cipher_mod, c->psi_cipher, sk2);

                        poly* term2 = new poly(res->ring_dim);
                        poly_handler::poly_mul_ntt(c->ciphertext[2], sk2, c->cipher_mod, c->psi_cipher, term2);

                        poly* term1 = new poly(res->ring_dim);
                        poly_handler::poly_mul_ntt(c->ciphertext[1], sk, c->cipher_mod, c->psi_cipher, term1);

                        poly_handler::poly_neg(term1, term1);
                        poly_handler::poly_add(c->ciphertext[0], term1, res);
                        poly_handler::poly_add(res, term2, res);
                        poly_handler::poly_mod(res, c->cipher_mod, res);
                        poly_handler::poly_bias_mod(res, c->cipher_mod, res);
                        poly_handler::poly_mod(res, c->plain_mod, res);

                        delete sk2;
                        delete term2;
                        delete term1;
                        break;
                    }
                    default:
                    {
                        return -1;
                    }
                }
                return 0;
            }
        }

        static bool eval_meta_eq(cipher* op1, cipher* op2)
        {
            if
            (
                op1->ring_dim != op2->ring_dim ||
                op1->level != op2->level ||
                op1->plain_mod != op2->plain_mod ||
                op1->cipher_mod != op2->cipher_mod ||
                op1->psi_plain != op2->psi_plain ||
                op1->psi_cipher != op2->psi_cipher
            )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        static int eval_add(cipher* op1, cipher* op2, cipher* res)
        {
            if((res->ring_dim != op1->ring_dim) || (res->ring_dim != op2->ring_dim))
            {
                return -1;
            }
            else if(!eval_meta_eq(op1, op2))
            {
                return -1;
            }
            else
            {
                switch(op1->level)
                {
                    case 1:
                    {
                        switch (res->level)
                        {
                            case 0:
                            {
                                res->level = 1;
                                poly* L1 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L1);
                                poly* L2 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L2);

                                break;
                            }
                            case 1:
                            {
                                break;
                            }
                            case 2:
                            {
                                res->ciphertext.erase(res->ciphertext.end() - 1);
                                res->level = 1;
                                
                                break;
                            }
                            default:
                            {
                                return -1;
                            }
                        }

                        for(int i = 0; i < 2; i++)
                        {
                            poly_handler::poly_add(op1->ciphertext[i], op2->ciphertext[i], res->ciphertext[i]);
                            poly_handler::poly_mod(res->ciphertext[i], res->cipher_mod, res->ciphertext[i]);
                        }
                        break;
                    }
                    case 2:
                    {
                        switch (res->level)
                        {
                            case 0:
                            {
                                res->level = 2;
                                poly* L1 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L1);
                                poly* L2 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L2);
                                poly* L3 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L3);
                                
                                
                                break;
                            }
                            case 1:
                            {
                                poly* L3 = new poly(res->ring_dim);
                                res->ciphertext.push_back(L3);
                                res->level = 2;

                                break;
                            }
                            case 2:
                            {
                                break;
                            }
                            default:
                            {
                                return -1;
                            }
                        }

                        for(int i = 0; i < 3; i++)
                        {
                            poly_handler::poly_add(op1->ciphertext[i], op2->ciphertext[i], res->ciphertext[i]);
                            poly_handler::poly_mod(res->ciphertext[i], res->cipher_mod, res->ciphertext[i]);
                        }
                                
                        break;
                    }
                    default:
                    {
                        return -1;
                    }
                }
                
                return 0;
            }
        }

        static int eval_mul(cipher* op1, cipher* op2, cipher* res)
        {
            if((res->ring_dim != op1->ring_dim) || (res->ring_dim != op2->ring_dim))
            {
                return -1;
            }
            else if(!eval_meta_eq(op1, op2))
            {
                return -1;
            }
            else if(op1->level != 1 || op2->level != 1)
            {
                return -1;
            }
            else
            {
                switch (res->level)
                {
                    case 0:
                    {
                        res->level = 2;
                        poly* L1 = new poly(res->ring_dim);
                        res->ciphertext.push_back(L1);
                        poly* L2 = new poly(res->ring_dim);
                        res->ciphertext.push_back(L2);
                        poly* L3 = new poly(res->ring_dim);
                        res->ciphertext.push_back(L3);

                        break;
                    }
                    case 1:
                    {
                        res->level = 2;
                        poly* L3 = new poly(res->ring_dim);
                        res->ciphertext.push_back(L3);

                        break;
                    }
                    case 2:
                    {
                        break;
                    }
                    default:
                    {
                        return -1;
                    }
                }

                cipher* clone1 = op1->clone();
                cipher* clone2 = op2->clone();

                poly_handler::poly_mul_ntt(clone1->ciphertext[0], clone2->ciphertext[0], res->cipher_mod, res->psi_cipher, res->ciphertext[0]);
                poly_handler::poly_mod(res->ciphertext[0], res->cipher_mod, res->ciphertext[0]);

                poly* temp = new poly(res->ring_dim);

                poly_handler::poly_mul_ntt(clone1->ciphertext[0], clone2->ciphertext[1], res->cipher_mod, res->psi_cipher, temp);
                poly_handler::poly_mul_ntt(clone1->ciphertext[1], clone2->ciphertext[0], res->cipher_mod, res->psi_cipher, res->ciphertext[1]);
                poly_handler::poly_add(temp, res->ciphertext[1], res->ciphertext[1]);
                poly_handler::poly_mod(res->ciphertext[1], res->cipher_mod, res->ciphertext[1]);
                poly_handler::poly_mul_ntt(clone1->ciphertext[1], clone2->ciphertext[1], res->cipher_mod, res->psi_cipher, res->ciphertext[2]);
                poly_handler::poly_mod(res->ciphertext[2], res->cipher_mod, res->ciphertext[2]);

                delete clone1;
                delete clone2;
                delete temp;
                
                return 0;
            }
        }
};

#endif