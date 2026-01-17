#ifndef CONTROL_H
#define CONTROL_H

#include "Struct.h"
#include "RLWE.h"

// ================================================================================ //
//                                      Plant                                       //
// ================================================================================ //

class plant
{
    public:
        std::vector<double*> A;
        std::vector<double*> B;
        std::vector<double*> C;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> u;
        int dim = 4;

        plant(std::vector<double> x_init)
        {
            double* A_row_1 = new double[4]{1, 0.0497730983685617, 0.00335251988390142, 0.00005577253479364};
            double* A_row_2 = new double[4]{0, 0.990924994598992, 0.134769006497267, 0.00335251988390142};
            double* A_row_3 = new double[4]{0, -0.000570156442840377, 1.03920568620351, 0.0506518405182903};
            double* A_row_4 = new double[4]{0, -0.0229198990641610, 1.57789260894128, 1.03920568620351};

            this->A.push_back(A_row_1);
            this->A.push_back(A_row_2);
            this->A.push_back(A_row_3);
            this->A.push_back(A_row_4);

            double* B_row_1 = new double[1]{0.00226901631438307};
            double* B_row_2 = new double[1]{0.0907500540100833};
            double* B_row_3 = new double[1]{0.00570156442840377};
            double* B_row_4 = new double[1]{0.229198990641610};

            this->B.push_back(B_row_1);
            this->B.push_back(B_row_2);
            this->B.push_back(B_row_3);
            this->B.push_back(B_row_4);

            double* C_row_1 = new double[4]{1, 0, 0, 0};
            double* C_row_2 = new double[4]{0, 0, 1, 0};

            this->C.push_back(C_row_1);
            this->C.push_back(C_row_2);

            this->x.resize(4);
            if(x_init.size() != 4)
            {
                for(int i = 0; i < 4; i++)
                {
                    this->x[i] = 0;
                }
            }
            else
            {
                this->x = x_init;
            }
            
            this->y.resize(2);
            this->u.resize(1);
        }

        ~plant()
        {
            for(auto& factor : this->A)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->A.clear();

            for(auto& factor : this->B)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->B.clear();

            for(auto& factor : this->C)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->C.clear();
        }

        void state_update(std::vector<double> u)
        {  
            std::vector<double> temp(4, 0);
            
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    temp[i] += this->A[i][j] * this->x[j];
                }
            }

            for(int i = 0; i < 4; i++)
            {
                temp[i] += this->B[i][0] * u[0];
            }

            this->x = temp;
        }

        void output()
        {
            this->y[0] = this->x[0];
            this->y[1] = this->x[2];
        }
};

// ================================================================================ //
//                                  Encrypted Controller                            //
// ================================================================================ //

class arx
{
    public:
        std::vector<cipher*> P_y;
        std::vector<cipher*> Q_u;
        std::vector<cipher*> mem_y;
        std::vector<cipher*> mem_u;
        std::vector<cipher*> temp1;
        std::vector<cipher*> temp2;
        cipher* temp;
        cipher* calc_res;

        arx()
        {
            this->P_y.resize(4);
            this->Q_u.resize(4);
            this->mem_y.resize(4);
            this->mem_u.resize(4);
            this->temp1.resize(4);
            this->temp2.resize(4);
        }

        ~arx()
        {
            for(auto& factor : this->P_y)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->P_y.clear();

            for(auto& factor : this->Q_u)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->Q_u.clear();

            for(auto& factor : this->mem_y)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_y.clear();

            for(auto& factor : this->mem_u)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->mem_u.clear();

            for(auto& factor : this->temp1)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->temp1.clear();

            for(auto& factor : this->temp2)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->temp2.clear();

            delete this->temp;
            delete this->calc_res;
        }

        void zero_set()
        {
            for(int i = 0; i < 4; i++)
            {
                this->mem_y[i] = this->temp->clone();
            }

            for(int i = 0; i < 4; i++)
            {
                this->mem_u[i] = this->temp->clone();
            }

            for(int i = 0; i < 4; i++)
            {
                this->temp1[i] = this->temp->clone();
            }

            for(int i = 0; i < 4; i++)
            {
                this->temp2[i] = this->temp->clone();
            }


            this->calc_res = this->temp->clone();
        }

        void calc()
        {
            for(int i = 0; i < 4; i++)
            {
                crypto_handler::eval_mul(this->P_y[i], this->mem_y[i], this->temp1[i]);
                crypto_handler::eval_mul(this->Q_u[i], this->mem_u[i], this->temp2[i]);
            }

            crypto_handler::eval_add(this->temp1[0], this->temp2[0], this->calc_res);
            for(int i = 1; i < 4; i++)
            {
                crypto_handler::eval_add(this->temp1[i], this->temp2[i], this->temp);
                crypto_handler::eval_add(this->calc_res, this->temp, this->calc_res);
            }
        }

        void mem_update(cipher* new_y, cipher* new_u)
        {
            for(int i = 0; i < 3; i++)
            {
                std::swap(this->mem_y[i], this->mem_y[i+1]);
                std::swap(this->mem_u[i], this->mem_u[i+1]);
            }
            this->mem_y[3]->ciphertext[0]->coeff = new_y->ciphertext[0]->coeff;
            this->mem_y[3]->ciphertext[1]->coeff = new_y->ciphertext[1]->coeff;
            this->mem_u[3]->ciphertext[0]->coeff = new_u->ciphertext[0]->coeff;
            this->mem_u[3]->ciphertext[1]->coeff = new_u->ciphertext[1]->coeff;
        }
};

// ================================================================================ //
//                                Verifiable Computation                            //
// ================================================================================ //

class vc
{
    public:
        std::vector<mr_cipher*> P_y;
        std::vector<mr_cipher*> Q_u;
        std::vector<cipher*> *mem_y;
        std::vector<cipher*> *mem_u;

        poly* r0;
        poly* r1;
        poly* s;
        gvec* g_r0;
        gvec* g_r1;
        gvec* g_s;

        poly* rz_f;
        poly* ro_f;
        poly* s_h;
        gvec* g_rz_f;
        gvec* g_ro_f;
        gvec* g_s_h;

        poly* rz_f_m_ro;
        poly* ro_f_m_rz;
        poly* s_h_m_rz;
        poly* s_h_m_ro;
        gvec* g_rz_f_m_ro;
        gvec* g_ro_f_m_rz;
        gvec* g_s_h_m_rz;
        gvec* g_s_h_m_ro;

        vc()
        {
            this->P_y.resize(4);
            this->Q_u.resize(4);
        }
        
        ~vc()
        {
            for(auto& factor : this->P_y)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->P_y.clear();

            for(auto& factor : this->Q_u)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->Q_u.clear();

            delete r0;
            delete r1;
            delete s;
            delete rz_f;
            delete ro_f;
            delete s_h;
            delete rz_f_m_ro;
            delete ro_f_m_rz;
            delete s_h_m_rz;
            delete s_h_m_ro;

            delete g_r0;
            delete g_r1;
            delete g_s;
            delete g_rz_f;
            delete g_ro_f;
            delete g_s_h;
            delete g_rz_f_m_ro;
            delete g_ro_f_m_rz;
            delete g_s_h_m_rz;
            delete g_s_h_m_ro;
        }

        void arx_coe_set(arx* ctrl)
        {
            for(int i = 0; i < 4; i++)
            {
                mr_cipher* PQ = new mr_cipher(ctrl->temp->ring_dim, ctrl->temp->plain_mod, ctrl->temp->cipher_mod);
                this->P_y[i] = PQ;
                this->Q_u[i] = PQ->clone();
                crypto_handler::cipher_2_mr_cipher(ctrl->P_y[i], this->P_y[i]);
                crypto_handler::cipher_2_mr_cipher(ctrl->Q_u[i], this->Q_u[i]);
            }
        }

        void link_memory(arx* ctrl)
        {
            this->mem_y = &ctrl->mem_y;
            this->mem_u = &ctrl->mem_u;
        }

        void set_ekf(poly* r0, poly* r1, poly* s)
        {
            this->r0 = r0;
            this->r1 = r1;
            this->s = s;

            this->rz_f = new poly(r0->ring_dim);
            this->ro_f = new poly(r1->ring_dim);
            this->s_h = new poly(r0->ring_dim);

            int N = r0->ring_dim / 8;

            for(int i = 0; i < 8; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    if(i == 0)
                    {
                        this->rz_f->coeff[j] = 0;
                        this->ro_f->coeff[j] = 0;
                    }
                    else if(i < 4)
                    {
                        this->rz_f->coeff[j + N * i] = this->r0->coeff[j + N * (i - 1)];
                        this->ro_f->coeff[j + N * i] = this->r1->coeff[j + N * (i - 1)];
                    }
                    else if(i == 4)
                    {
                        this->rz_f->coeff[j + N * i] = 0;
                        this->ro_f->coeff[j + N * i] = 0;
                    }
                    else
                    {
                        this->rz_f->coeff[j + N * i] = this->r0->coeff[j + N * (i - 1)];
                        this->ro_f->coeff[j + N * i] = this->r1->coeff[j + N * (i - 1)];
                    }
                }   
            }

            poly_handler::poly_set_2_zero(this->s_h);
            for(int i = 0; i < 8; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    if(j < (N / 3))
                    {
                        for(int k = 0; k < (N / 3); k++)
                        {
                            if(i < 4)
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->P_y[i]->ciphertext[0]->entry[j + this->P_y[i]->ring_dim * k];
                            }
                            else
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->Q_u[i - 4]->ciphertext[0]->entry[j + this->Q_u[i - 4]->ring_dim * k];
                            }
                            
                        }
                        for(int k = (N / 3); k < 2 * (N / 3); k++)
                        {
                            if(i < 4)
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->P_y[i]->ciphertext[1]->entry[j + this->P_y[i]->ring_dim * (k - (N / 3))];
                            }
                            else
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->Q_u[i - 4]->ciphertext[1]->entry[j + this->Q_u[i - 4]->ring_dim * (k - (N / 3))];
                            }
                        }
                    }
                    else if(j >= (N / 3) && j < 2 * (N / 3))
                    {
                        for(int k = (N / 3); k < 2 * (N / 3); k++)
                        {
                            if(i < 4)
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->P_y[i]->ciphertext[0]->entry[(j - (N / 3)) + this->P_y[i]->ring_dim * (k - (N / 3))];
                            }
                            else
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->Q_u[i - 4]->ciphertext[0]->entry[(j - (N / 3)) + this->Q_u[i - 4]->ring_dim * (k - (N / 3))];
                            }
                        }
                    }
                    else
                    {
                        for(int k = 2 * (N / 3); k < 3 * (N / 3); k++)
                        {
                            if(i < 4)
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->P_y[i]->ciphertext[1]->entry[(j - 2 * (N / 3)) + this->P_y[i]->ring_dim * (k -  2 * (N / 3))];
                            }
                            else
                            {
                                this->s_h->coeff[N * i + j] += this->s->coeff[k] * this->Q_u[i - 4]->ciphertext[1]->entry[(j - 2 * (N / 3)) + this->Q_u[i - 4]->ring_dim * (k -  2 * (N / 3))];
                            }
                        }
                    }
                }
            }
            
            poly* neg_rz = this->r0->clone();
            poly* neg_ro = this->r1->clone();

            this->rz_f_m_ro = new poly(r0->ring_dim);
            this->ro_f_m_rz = new poly(r1->ring_dim);
            this->s_h_m_rz = new poly(s->ring_dim);
            this->s_h_m_ro = new poly(s->ring_dim);

            poly_handler::poly_add(this->rz_f, neg_ro, this->rz_f_m_ro);
            poly_handler::poly_add(this->ro_f, neg_rz, this->ro_f_m_rz);
            poly_handler::poly_add(this->s_h, neg_rz, this->s_h_m_ro);
            poly_handler::poly_add(this->s_h, neg_ro, this->s_h_m_rz);

            for(auto& factor : this->P_y)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->P_y.clear();

            for(auto& factor : this->Q_u)
            {
                if(factor != nullptr)
                {
                    delete factor;
                    factor = nullptr;
                }
            }
            this->Q_u.clear();

            delete neg_rz;
            delete neg_ro;
        }

        void up_2_g(mpz_class g_mod, mpz_class g_gen)
        {
            this->g_r0 = new gvec(this->r0->ring_dim, g_mod, g_gen);
            this->g_r1 = new gvec(this->r1->ring_dim, g_mod, g_gen);
            this->g_s = new gvec(this->s->ring_dim, g_mod, g_gen);

            group_handler::poly_2_gvec(this->r0, this->g_r0);
            group_handler::poly_2_gvec(this->r1, this->g_r1);
            group_handler::poly_2_gvec(this->s, this->g_s);

            this->g_rz_f = new gvec(this->rz_f->ring_dim, g_mod, g_gen);
            this->g_ro_f = new gvec(this->ro_f->ring_dim, g_mod, g_gen);
            this->g_s_h = new gvec(this->s_h->ring_dim, g_mod, g_gen);

            group_handler::poly_2_gvec(this->rz_f, this->g_rz_f);
            group_handler::poly_2_gvec(this->ro_f, this->g_ro_f);
            group_handler::poly_2_gvec(this->s_h, this->g_s_h);

            this->g_rz_f_m_ro = new gvec(this->rz_f_m_ro->ring_dim, g_mod, g_gen);
            this->g_ro_f_m_rz = new gvec(this->ro_f_m_rz->ring_dim, g_mod, g_gen);
            this->g_s_h_m_rz = new gvec(this->s_h_m_rz->ring_dim, g_mod, g_gen);
            this->g_s_h_m_ro = new gvec(this->s_h_m_ro->ring_dim, g_mod, g_gen);

            group_handler::poly_2_gvec(this->rz_f_m_ro, this->g_rz_f_m_ro);
            group_handler::poly_2_gvec(this->ro_f_m_rz, this->g_ro_f_m_rz);
            group_handler::poly_2_gvec(this->s_h_m_rz, this->g_s_h_m_rz);
            group_handler::poly_2_gvec(this->s_h_m_ro, this->g_s_h_m_ro);
        }
};

#endif