#ifndef GROUP_H
#define GROUP_H

#include "Struct.h"

// ================================================================================ //
//                                  Group Structure                                 //
// ================================================================================ //

class gvec
{
    public:
        std::vector<mpz_class> vec;
        int size;
        mpz_class g_mod;
        mpz_class g_gen;

    gvec(int size, mpz_class g_mod, mpz_class g_gen):
        size(size),
        g_mod(g_mod),
        g_gen(g_gen)
    {
        this->vec.resize(size);
    }

    ~gvec()
    {
        this->vec.clear();
    }

    gvec* clone()
    {
        gvec* clone = new gvec(this->size, this->g_mod, this->g_gen);
        clone->vec = this->vec;
        
        return clone;
    }
};

// ================================================================================ //
//                                  Group Handler                                   //
// ================================================================================ //

class group_handler
{
    public:
        static int poly_2_gvec(poly* op, gvec* res)
        {
            if(op->ring_dim != res->size)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->size; i++)
                {
                    
                    mpz_powm(res->vec[i].get_mpz_t(), res->g_gen.get_mpz_t(), op->coeff[i].get_mpz_t() , res->g_mod.get_mpz_t());
                }

                return 0;
            }
        }

        static int gvec_2_poly(gvec* op, const mpz_class& mod, poly* res)
        {
            if(op->size != res->ring_dim)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    
                    res->coeff[i] = op->vec[i] % mod;
                }

                return 0;
            }
        }

        static bool eval_meta_eq(gvec* op1, gvec* op2)
        {
            if
            (
                op1->size != op2->size ||
                op1->g_mod != op2->g_mod ||
                op1->g_gen != op2->g_gen
            )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        static int group_exp(gvec* op1, poly* op2, gvec* res)
        {
            if((res->size != op1->size) || (res->size != op2->ring_dim))
            {
                return -1;
            }
            else
            {
                if(!eval_meta_eq(op1, res))
                {
                    return -1;
                }
                else
                {
                    for(int i = 0; i < res->size; i++)
                    {
                        mpz_powm(res->vec[i].get_mpz_t(), res->g_gen.get_mpz_t(), op2->coeff[i].get_mpz_t(), res->g_mod.get_mpz_t());
                    }

                    return 0;
                }
            }
        }

        static int group_mul(gvec* op1, gvec* op2, gvec* res)
        {
            if((res->size != op1->size) || (res->size != op2->size))
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->size; i++)
                {
                    res->vec[i] = (op1->vec[i] * op2->vec[i]) % res->g_mod;
                }

                return 0;
            }
        }

        static int group_dot(gvec* op1, gvec* op2, mpz_class& res)
        {
            if(op1->size != op2->size)
            {
                return -1;
            }
            else
            {
                if(op1->g_mod != op2->g_mod)
                {
                    return -1;
                }
                else
                {
                    res = 1;

                    for(int i = 0; i < op1->size; i++)
                    {
                        res = (res + (op1->vec[i] * op2->vec[i]) % op1->g_mod) % op1->g_mod;
                    }

                    return 0;
                }   
            }
        }
};

#endif