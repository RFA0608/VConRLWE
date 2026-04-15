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

struct jpoint
{
    mpz_class x;
    mpz_class y;
    mpz_class z;

    jpoint() : x(0), y(0), z(0) {}
};

class eccvec
{
    public:
        std::vector<point> p;
        std::vector<jpoint> jp;
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
        this->jp.resize(size);
    }

    ~eccvec()
    {
        this->p.clear();
        this->jp.clear();
    }

    eccvec* clone()
    {
        eccvec* clone = new eccvec(this->size);
        clone->p = this->p;
        clone->jp = this->jp;
        
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

                        res->jp[i].x = 0;
                        res->jp[i].y = 0;
                        res->jp[i].z = 0;
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
                        
                        res->jp[i].x = res->p[i].x;
                        res->jp[i].y = res->p[i].y;
                        res->jp[i].z = 1;
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

        static int j_double(jpoint& p, jpoint& r, mpz_class& ecc_mod)
        {
            if(p.z == 0)
            {
                r.z == 0;
            }
            else
            {
                mpz_class Y_sq = (p.y * p.y) % ecc_mod;
                mpz_class S = (4 * p.x * Y_sq) % ecc_mod;
                mpz_class M = (3 * p.x * p.x) % ecc_mod;
                
                mpz_class X3 = ((M * M - 2 * S) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class Y_4 = (Y_sq * Y_sq) % ecc_mod;
                
                mpz_class temp = ((S - X3) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class Y3 = ((M * temp - 8 * Y_4) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class Z3 = (2 * p.y * p.z) % ecc_mod;

                r.x = X3; 
                r.y = Y3; 
                r.z = Z3;
            }
            return 0;
        }

        static int j_add(jpoint& p1, jpoint& p2, jpoint& r, mpz_class& ecc_mod)
        {
            if(p1.z == 0)
            {
                r = p2;
            }
            else if(p2.z == 0)
            {
                r = p1;
            }
            else
            {
                mpz_class Z1_sq = (p1.z * p1.z) % ecc_mod;
                mpz_class Z2_sq = (p2.z * p2.z) % ecc_mod;
                mpz_class Z1_cu = (Z1_sq * p1.z) % ecc_mod;
                mpz_class Z2_cu = (Z2_sq * p2.z) % ecc_mod;

                mpz_class U1 = (p1.x * Z2_sq) % ecc_mod;
                mpz_class U2 = (p2.x * Z1_sq) % ecc_mod;
                mpz_class S1 = (p1.y * Z2_cu) % ecc_mod;
                mpz_class S2 = (p2.y * Z1_cu) % ecc_mod;

                if(U1 == U2)
                {
                    if(S1 == S2) 
                    { 
                        j_double(p1, r, ecc_mod); 
                        return 0; 
                    }
                    else 
                    { 
                        r.x = 0;
                        r.y = 0;
                        return 0; 
                    }
                }

                mpz_class H = ((U2 - U1) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class R = ((S2 - S1) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class H_sq = (H * H) % ecc_mod;
                mpz_class H_cu = (H * H_sq) % ecc_mod;
                mpz_class U1_H_sq = (U1 * H_sq) % ecc_mod;

                mpz_class X3 = ((R * R - H_cu - 2 * U1_H_sq) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class temp1 = ((U1_H_sq - X3) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class temp2 = (S1 * H_cu) % ecc_mod;
                mpz_class Y3 = ((R * temp1 - temp2) % ecc_mod + ecc_mod) % ecc_mod;
                mpz_class Z3 = (((H * p1.z) % ecc_mod) * p2.z) % ecc_mod;

                r.x = X3; 
                r.y = Y3; 
                r.z = Z3;
            }
            return 0;
        }

        static int to_affine(jpoint& p, point& r, mpz_class& ecc_mod)
        {
            if(p.z == 0)
            {
                r.x = 0; 
                r.y = 0;
                
            }
            else
            {
                 mpz_class Z_inv, Z_inv_sq, Z_inv_cu;
                mpz_invert(Z_inv.get_mpz_t(), p.z.get_mpz_t(), ecc_mod.get_mpz_t());
                
                Z_inv_sq = (Z_inv * Z_inv) % ecc_mod;
                Z_inv_cu = (Z_inv_sq * Z_inv) % ecc_mod;

                r.x = ((p.x * Z_inv_sq) % ecc_mod + ecc_mod) % ecc_mod;
                r.y = ((p.y * Z_inv_cu) % ecc_mod + ecc_mod) % ecc_mod;
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

        static int jpoint_mul(jpoint& p, mpz_class& n, jpoint& r, mpz_class& ecc_mod)
        {
            r.x = p.x;
            r.y = p.y;
            r.z = 1;

            jpoint D;

            for(int j = 0; j < 256; j++)
            {
                if(mpz_tstbit(n.get_mpz_t(), j))
                {
                    j_add(r, D, r, ecc_mod);
                }
                j_double(D, D, ecc_mod);
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

        static int ecc_dot_j(eccvec* P, poly* op, point& res)
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

                        R->jp[i].x = 0;
                        R->jp[i].y = 0;
                        R->jp[i].z = 0;
                    }
                    else if(P->p[i].x == 0 && P->p[i].y == 0)
                    {
                        R->p[i].x = 0;
                        R->p[i].y = 0;

                        R->jp[i].x = 0;
                        R->jp[i].y = 0;
                        R->jp[i].z = 0;
                    }
                    else
                    {
                        jpoint D;
                        D.x = P->jp[i].x;
                        D.y = P->jp[i].y;
                        D.z = P->jp[i].z;
                        for(int j = 0; j < 256; j++)
                        {
                            if(mpz_tstbit(op->coeff[i].get_mpz_t(), j))
                            {
                                j_add(R->jp[i], D, R->jp[i], R->ecc_mod);

                            }
                            j_double(D, D, R->ecc_mod);
                        }

                        to_affine(R->jp[i], R->p[i], R->ecc_mod);
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