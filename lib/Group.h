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
                #pragma omp parallel for
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

        static int group_exp(poly* op, gvec* res)
        {
            if(res->size != op->ring_dim)
            {
                return -1;
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i < res->size; i++)
                {
                    mpz_powm(res->vec[i].get_mpz_t(), res->g_gen.get_mpz_t(), op->coeff[i].get_mpz_t(), res->g_mod.get_mpz_t());
                }

                return 0;
            }
        }

        // static void group_scal_exp(mpz_class& op1, mpz_class& op2, mpz_class& g_mod, mpz_class& res)
        // {
        //     mpz_powm(res.get_mpz_t(), op1.get_mpz_t(), op2.get_mpz_t(), g_mod.get_mpz_t());
        // }

        static int group_mul(gvec* op1, gvec* op2, gvec* res)
        {
            if(!group_handler::eval_meta_eq(op1, res) || !group_handler::eval_meta_eq(op2, res))
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

        // static int group_dot(gvec* op1, poly* op2, mpz_class& res)
        // {
        //     if(op1->size != op2->ring_dim)
        //     {
        //         return -1;
        //     }
        //     else
        //     {
        //         std::vector<mpz_class> temp(op1->size);
        //         res = 1;

        //         #pragma omp parallel for
        //         for(int i = 0; i < op1->size; i++)
        //         {
        //             if(op2->coeff[i] == 0)
        //             {
        //                 temp[i] = 1;
        //             }
        //             else if(op2->coeff[i] == 1)
        //             {
        //                 temp[i] = op1->vec[i] %  op1->g_mod;
        //             }
        //             else
        //             {
        //                 mpz_powm(temp[i].get_mpz_t(), op1->vec[i].get_mpz_t(), op2->coeff[i].get_mpz_t(), op1->g_mod.get_mpz_t());
        //             }   
        //         }

        //         for(int i = 0; i < op1->size; i++)
        //         {
        //             res *= temp[i];
        //             res %= op1->g_mod;
        //         }

        //         return 0;
        //     }
        // }

        static int group_dot(gvec* op1, poly* op2, mpz_class& res)
        {
            if(op1->size != op2->ring_dim) 
            {
                return -1;
            }
            
            res = 1;
            int num_threads = omp_get_max_threads();
            std::vector<mpz_class> thread_results(num_threads, 1);

            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                mpz_class local_res = 1;
                mpz_class p_res;

                #pragma omp for
                for(int i = 0; i < op1->size; i++)
                {
                    const mpz_class& exp = op2->coeff[i];
                    
                    if(exp == 0) continue;
                    
                    if(exp == 1) {
                        local_res = (local_res * op1->vec[i]) % op1->g_mod;
                    }
                    else {
                        mpz_powm(p_res.get_mpz_t(), op1->vec[i].get_mpz_t(), exp.get_mpz_t(), op1->g_mod.get_mpz_t());
                        local_res = (local_res * p_res) % op1->g_mod;
                    }
                }
                thread_results[tid] = local_res;
            }

            for(int i = 0; i < num_threads; i++) {
                res = (res * thread_results[i]) % op1->g_mod;
            }

            return 0;
        }
};

#endif
