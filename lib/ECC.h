#ifndef ECC_H
#define ECC_H

#include "Struct.h"
#include <omp.h>

// ================================================================================ //
//                                  Group Structure                                 //
// ================================================================================ //

class eccvec
{
    public:
        std::vector<mpz_class> x;
        std::vector<mpz_class> y;
        int size;
        mpz_class ecc_mod; 
        mpz_class Gx;
        mpz_class Gy;
        
        // y^2 = x^3 + 5 (0, 0 is infty origin)
        mpz_class a = 0;
        mpz_class b = 4; 
        
    eccvec(int size):
        size(size)
    {
        this->ecc_mod = "28948022309329048855892746252171976963363056481941560715954676764349967630337";
        this->Gx = "28948022309329048855892746252171976963363056481941560715954676764349967630336";
        this->Gy = "2";
        this->x.resize(size);
        this->y.resize(size);
    }

    ~eccvec()
    {
        this->x.clear();
        this->y.clear();
    }

    eccvec* clone()
    {
        eccvec* clone = new eccvec(this->size);
        clone->x = this->x;
        clone->y = this->y;
        
        return clone;
    }
};

// ================================================================================ //
//                                  Group Handler                                   //
// ================================================================================ //

class ecc_handler
{
    public:
        static int poly_2_ecc(poly* op, eccvec* table, eccvec* res)
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
                    if(op->coeff[i] == 0)
                    {
                        res->x[i] = 0;
                        res->y[i] = 0;
                    }
                    else
                    {
                        res->x[i] = 0;
                        res->y[i] = 0;
                        for(int j = 0; j < 256; j++)
                        {
                            
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                point_add(res->x[i], res->y[i], table->x[j], table->y[j], res->x[i], res->y[i], res->ecc_mod);
                            }
                        }
                    }
                } 

                return 0;
            }
        };

        static int table_set(eccvec* table)
        {
            for(int i = 0; i < 256; i++)
            {
                if(i == 0)
                {
                    table->x[i] = table->Gx;
                    table->y[i] = table->Gy;
                }
                else
                {
                    mpz_class lambda;
                    mpz_class temp;
                    temp = ((2 * table->y[i-1]) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                    mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), table->ecc_mod.get_mpz_t());
                    lambda = lambda * (3 * table->x[i-1] * table->x[i-1] + table->a);

                    table->x[i] = ((lambda * lambda - table->x[i-1] - table->x[i-1]) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                    table->y[i] = ((lambda * (table->x[i-1] - table->x[i]) - table->y[i-1]) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                }
            }
            return 0;
        }

        static int ecc_add(eccvec* P, eccvec* Q, eccvec* R)
        {
            if((P->size != R->size) || (Q->size != R->size))
            {
                return -1;
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i < R->size; i++)
                {
                    mpz_class lambda;
                    if((P->x[i] == 0 && P->y[i]) == 0 && (Q->x[i] == 0 && Q->y[i] == 0))
                    {
                        R->x[i] = 0;
                        R->y[i] = 0;
                    }
                    else if(P->x[i] == 0 && P->y[i] == 0)
                    {
                        R->x[i] = Q->x[i];
                        R->y[i] = Q->y[i];
                    }
                    else if(Q->x[i] == 0 && Q->y[i] == 0)
                    {
                        R->x[i] = P->x[i];
                        R->y[i] = P->y[i];
                    }
                    else if(P->x[i] == Q->x[i])
                    {
                        R->x[i] = 0;
                        R->y[i] = 0;
                    }
                    else if(P->x[i] != Q->x[i] || P->y[i] != Q->y[i])
                    {
                        mpz_class temp;
                        temp = ((Q->x[i] - P->x[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                        mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), R->ecc_mod.get_mpz_t());
                        lambda = lambda * (Q->y[i] - P->y[i]);

                        R->x[i] = ((lambda * lambda - P->x[i] - Q->x[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                        R->y[i] = ((lambda * (P->x[i] - R->x[i]) - P->y[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                    }
                    else
                    {
                        mpz_class temp;
                        temp = ((2 * P->y[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                        mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), R->ecc_mod.get_mpz_t());
                        lambda = lambda * (3 * P->x[i] * P->x[i] + R->a);

                        R->x[i] = ((lambda * lambda - P->x[i] - Q->x[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                        R->y[i] = ((lambda * (P->x[i] - R->x[i]) - P->y[i]) % R->ecc_mod + R->ecc_mod) % R->ecc_mod;
                    }
                }

                return 0;
            }
        };

        static int point_add(mpz_class& x1, mpz_class& y1, mpz_class& x2, mpz_class& y2, mpz_class& rx, mpz_class& ry, mpz_class& ecc_mod)
        {
            if((x1 == 0 && y1 == 0) && (x2 == 0 && y2 == 0))
            {
                rx = 0;
                ry = 0;
            }
            else if(x1 == 0 && y1 == 0)
            {
                rx = x2;
                ry = y2;
            }
            else if(x2 == 0 && y2 == 0)
            {
                rx = x1;
                ry = y2;
            }
            else if(x1 == x2 && y1 != y2)
            {
                rx = 0;
                ry = 0;
            }
            else if(x1 != x2 || y1 != y2)
            {
                mpz_class lambda;
                mpz_class temp, x1c, y1c, x2c, y2c;
                x1c = x1; y1c = y1;
                x2c = x2; y2c = y2;
                temp = ((x2c - x1c) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), ecc_mod.get_mpz_t());
                lambda = lambda * (y2c - y1c);

                rx = ((lambda * lambda - x1c - x2c) % ecc_mod + ecc_mod) % ecc_mod;
                ry = ((lambda * (x1c - rx) - y1c) % ecc_mod + ecc_mod) % ecc_mod;
            }
            else
            {
                mpz_class lambda;
                mpz_class temp, x1c, y1c, x2c, y2c;
                x1c = x1; y1c = y1;
                x2c = x2; y2c = y2;
                temp = ((2 * y1c) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(),ecc_mod.get_mpz_t());
                lambda = lambda * 3 * (x1c * x1c); // (x1c * x1c + a) but a = 0 in this struct

                rx = ((lambda * lambda - x1c - x2c) % ecc_mod + ecc_mod) % ecc_mod;
                ry = ((lambda * (x1c - rx) - y1c) % ecc_mod + ecc_mod) % ecc_mod;
            }

            return 0;
        }

        static int ecc_mul(eccvec* P, poly* op, eccvec* R)
        {
            if((P->size != R->size) || (op->ring_dim != R->size))
            {
                return -1;
            }
            else
            {
                #pragma omp parallel for
                for(int i = 0; i < R->size; i++)
                {
                    if(op->coeff[i] == 0)
                    {
                        R->x[i] = 0;
                        R->y[i] = 0;
                    }
                    else if(P->x[i] == 0 && P->y[i] == 0)
                    {
                        R->x[i] = 0;
                        R->y[i] = 0;
                    }
                    else
                    {
                        mpz_class Dx = P->x[i];
                        mpz_class Dy = P->y[i];
                        R->x[i] = 0;
                        R->y[i] = 0;
                        for(int j = 0; j < 256; j++)
                        {
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                point_add(R->x[i], R->y[i], Dx, Dy, R->x[i], R->y[i], R->ecc_mod);

                            }
                            point_add(Dx, Dy, Dx, Dy, Dx, Dy, R->ecc_mod);
                        }
                    }
                } 
                
                return 0;
            }
        }
};

#endif