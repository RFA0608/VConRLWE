#ifndef ECC_H
#define ECC_H

#include "Struct.h"
#include <omp.h>

// ================================================================================ //
//                                   ECC Structure                                  //
// ================================================================================ //

struct point
{
    mpz_class x;
    mpz_class y;

    point() : x(0), y(0) {}
};

class eccvec
{
    public:
        std::vector<point> p;
        int size;
        mpz_class ecc_mod; 
        point G;
        
        // y^2 = x^3 + 5 (0, 0 is infty origin)
        mpz_class a = 0;
        mpz_class b = 5; 
        
    eccvec(int size):
        size(size)
    {
        this->ecc_mod = "28948022309329048855892746252171976963363056481941560715954676764349967630337";
        this->G.x = "28948022309329048855892746252171976963363056481941560715954676764349967630336";
        this->G.y = "2";
        this->p.resize(size);
    }

    ~eccvec()
    {
        this->p.clear();
    }

    eccvec* clone()
    {
        eccvec* clone = new eccvec(this->size);
        clone->p = this->p;
        
        return clone;
    }
};

// ================================================================================ //
//                                   ECC Handler                                    //
// ================================================================================ //

class ecc_handler
{
    public:
        static bool compare(point& p, point& q)
        {
            if((p.x == q.x) && (p.y == q.y))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

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
                        res->p[i].x = 0;
                        res->p[i].y = 0;
                    }
                    else
                    {
                        res->p[i].x = 0;
                        res->p[i].y = 0;
                        for(int j = 0; j < 256; j++)
                        {
                            
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                point_add(res->p[i], table->p[j], res->p[i], res->ecc_mod);
                            }
                        }
                    }
                } 

                return 0;
            }
        };

        static int table_set(eccvec* table)
        {
            #pragma omp parallel for
            for(int i = 0; i < 256; i++)
            {
                if(i == 0)
                {
                    table->p[i].x = table->G.x;
                    table->p[i].y = table->G.y;
                }
                else
                {
                    mpz_class lambda;
                    mpz_class temp;
                    temp = ((2 * table->p[i-1].y) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                    mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), table->ecc_mod.get_mpz_t());
                    // lambda = lambda * (3 * table->p[i-1].x * table->p[i-1].x + table->a);
                    lambda = ((lambda * (3 * table->p[i-1].x * table->p[i-1].x + table->a)) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;

                    table->p[i].x = ((lambda * lambda - table->p[i-1].x - table->p[i-1].x) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                    table->p[i].y = ((lambda * (table->p[i-1].x - table->p[i].x) - table->p[i-1].y) % table->ecc_mod + table->ecc_mod) % table->ecc_mod;
                }
            }
            return 0;
        }

        static int point_add(point& p1, point& p2, point& r, mpz_class& ecc_mod)
        {
            if((p1.x == 0 && p1.y == 0) && (p2.x == 0 && p2.y == 0))
            {
                r.x = 0;
                r.y = 0;
            }
            else if(p1.x == 0 && p1.y == 0)
            {
                r.x = p2.x;
                r.y = p2.y;
            }
            else if(p2.x == 0 && p2.y == 0)
            {
                r.x = p1.x;
                r.y = p1.y;
            }
            else if(p1.x == p2.x && p1.y != p2.y)
            {
                r.x = 0;
                r.y = 0;
            }
            else if(p1.x != p2.x || p1.y != p2.y)
            {
                mpz_class lambda;
                mpz_class temp, x1c, y1c, x2c, y2c;
                x1c = p1.x; y1c = p1.y;
                x2c = p2.x; y2c = p2.y;
                temp = ((x2c - x1c) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(), ecc_mod.get_mpz_t());
                // lambda = lambda * (y2c - y1c);
                lambda = ((lambda * (y2c - y1c)) % ecc_mod + ecc_mod) % ecc_mod;

                r.x = ((lambda * lambda - x1c - x2c) % ecc_mod + ecc_mod) % ecc_mod;
                r.y = ((lambda * (x1c - r.x) - y1c) % ecc_mod + ecc_mod) % ecc_mod;
            }
            else
            {
                mpz_class lambda;
                mpz_class temp, x1c, y1c, x2c, y2c;
                x1c = p1.x; y1c = p1.y;
                x2c = p2.x; y2c = p2.y;
                temp = ((2 * y1c) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_invert(lambda.get_mpz_t(), temp.get_mpz_t(),ecc_mod.get_mpz_t());
                // lambda = lambda * 3 * (x1c * x1c); // (x1c * x1c + a) but a = 0 in this struct
                lambda = ((lambda * 3 * (x1c * x1c)) % ecc_mod + ecc_mod) % ecc_mod;

                r.x = ((lambda * lambda - x1c - x2c) % ecc_mod + ecc_mod) % ecc_mod;
                r.y = ((lambda * (x1c - r.x) - y1c) % ecc_mod + ecc_mod) % ecc_mod;
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
                    point_add(P->p[i], Q->p[i], R->p[i], R->ecc_mod);
                }

                return 0;
            }
        }

        static int point_mul(point& p, mpz_class& n, point& r, mpz_class& ecc_mod)
        { 
            r.x = 0;
            r.y = 0;
            point D;
            D.x = p.x;
            D.y = p.y;
            for(int j = 0; j < 256; j++)
            {
                if(mpz_tstbit(n.get_mpz_t(), j))
                {
                    point_add(r, D, r, ecc_mod);

                }
                point_add(D, D, D, ecc_mod);
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
                        R->p[i].x = 0;
                        R->p[i].y = 0;
                    }
                    else if(P->p[i].x == 0 && P->p[i].y == 0)
                    {
                        R->p[i].x = 0;
                        R->p[i].y = 0;
                    }
                    else
                    {
                        R->p[i].x = 0;
                        R->p[i].y = 0;
                        point D;
                        D.x = P->p[i].x;
                        D.y = P->p[i].y;
                        for(int j = 0; j < 256; j++)
                        {
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                point_add(R->p[i], D, R->p[i], R->ecc_mod);

                            }
                            point_add(D, D, D, R->ecc_mod);
                        }
                    }
                } 
                
                return 0;
            }
        }

        static int ecc_dot(eccvec* P, poly* op, point& res)
        {
            if(P->size != op->ring_dim)
            {
                return -1;
            }
            else
            {
                eccvec* R = new eccvec(P->size);

                #pragma omp parallel for
                for(int i = 0; i < R->size; i++)
                {
                    if(op->coeff[i] == 0)
                    {
                        R->p[i].x = 0;
                        R->p[i].y = 0;
                    }
                    else if(P->p[i].x == 0 && P->p[i].y == 0)
                    {
                        R->p[i].x = 0;
                        R->p[i].y = 0;
                    }
                    else
                    {
                        point D;
                        D.x = P->p[i].x;
                        D.y = P->p[i].y;
                        for(int j = 0; j < 256; j++)
                        {
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                point_add(R->p[i], D, R->p[i], R->ecc_mod);

                            }
                            point_add(D, D, D, R->ecc_mod);
                        }
                    }
                } 
                
                res.x = 0;
                res.y = 0;
                for(int i = 0; i < R->size; i++)
                {
                    point_add(res, R->p[i], res, R->ecc_mod);
                }
                
                delete R;
                return 0;
            }
        }
};

#endif