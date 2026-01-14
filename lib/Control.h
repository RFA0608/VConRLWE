#ifndef CONTROL_H
#define CONTROL_H

#include "Struct.h"
#include "RLWE.h"

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

            this->x.resize(this->dim);
            if(x_init.size() != this->dim)
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
        int dim = 4;

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
            for(int i = 0; i < this->dim; i++)
            {
                this->mem_y[i] = this->temp->clone();
            }

            for(int i = 0; i < this->dim; i++)
            {
                this->mem_u[i] = this->temp->clone();
            }

            for(int i = 0; i < this->dim; i++)
            {
                this->temp1[i] = this->temp->clone();
            }

            for(int i = 0; i < this->dim; i++)
            {
                this->temp2[i] = this->temp->clone();
            }


            this->calc_res = this->temp->clone();
        }

        void calc()
        {
            for(int i = 0; i < this->dim; i++)
            {
                crypto_handler::eval_mul(this->P_y[i], this->mem_y[i], this->temp1[i]);
                crypto_handler::eval_mul(this->Q_u[i], this->mem_u[i], this->temp2[i]);
            }

            crypto_handler::eval_add(this->temp1[0], this->temp2[0], this->calc_res);
            for(int i = 1; i < this->dim; i++)
            {
                crypto_handler::eval_add(this->temp1[i], this->temp2[i], this->temp);
                crypto_handler::eval_add(this->calc_res, this->temp, this->calc_res);
            }
        }

        void mem_update(cipher* new_y, cipher* new_u)
        {
            for(int i = 0; i < this->dim-1; i++)
            {
                std::swap(this->mem_y[i], this->mem_y[i+1]);
                std::swap(this->mem_u[i], this->mem_u[i+1]);
            }
            this->mem_y[dim-1]->ciphertext[0]->coeff = new_y->ciphertext[0]->coeff;
            this->mem_y[dim-1]->ciphertext[1]->coeff = new_y->ciphertext[1]->coeff;
            this->mem_u[dim-1]->ciphertext[0]->coeff = new_u->ciphertext[0]->coeff;
            this->mem_u[dim-1]->ciphertext[1]->coeff = new_u->ciphertext[1]->coeff;
        }
};


#endif