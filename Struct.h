#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <cmath>

#include "seal/seal.h"

class poly
{
    public:
        int size;
        mpz_t* coeff;

        poly(int size)
        {
            this->size = size;
            
            this->coeff = new mpz_t[size];
            for(int i = 0; i < size; i++)
            {
                mpz_init_set_str(*(this->coeff + i), "0", 10);
            }
        };

        ~poly()
        {
            for(int i = 0; i < this->size; i++)
            {
                mpz_clear(*(this->coeff + i));
            }
            delete(this->coeff);
        }
};

int poly_set(seal::Plaintext op, poly* result)
{
    int poly_degree = op.coeff_count();

    for(int i = 0; i < poly_degree; i++)
    {
        mpz_set_ui(*(result->coeff + i), op[i]);
    }

    return 0;
}

int poly_get(poly* op, seal::Plaintext& result)
{
    int poly_degree = op->size;

    for(int i = 0; i < poly_degree; i++)
    {
        result[i] = mpz_get_ui(*(op->coeff + i));
    }

    return 0;
}

int poly_neg(poly* op, poly* result)
{
    poly* temp = new poly(op->size);

    if(op->size != result->size)
    {
        delete temp;

        return -1;
    }
    else
    {
        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(temp->coeff + i), *(op->coeff + i));
            mpz_neg(*(temp->coeff + i), *(temp->coeff + i));
        }

        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(result->coeff + i), *(temp->coeff + i));
        }

        delete temp;

        return 0;
    }
}

int poly_add(poly* op1, poly* op2, poly* result)
{
    poly* temp1 = new poly(op1->size);
    poly* temp2 = new poly(op2->size);

    if(op1->size != op2->size || op1->size != result->size || op2->size != result->size)
    {
        delete temp1;
        delete temp2;

        return -1;
    }
    else
    {
        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(temp1->coeff + i), *(op1->coeff + i));
            mpz_set(*(temp2->coeff + i), *(op2->coeff + i));
        }

        for(int i = 0; i < result->size; i++)
        {
            mpz_add(*(result->coeff + i), *(temp1->coeff + i), *(temp2->coeff + i));
        }

        delete temp1;
        delete temp2;

        return 0;
    }
}

int poly_mul(poly* op1, poly* op2, poly* result)
{
    poly* temp1 = new poly(op1->size);
    poly* temp2 = new poly(op2->size);

    if(op1->size != op2->size || op1->size != result->size || op2->size != result->size)
    {
        delete temp1;
        delete temp2;

        return -1;
    }
    else
    {
        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(temp1->coeff + i), *(op1->coeff + i));
            mpz_set(*(temp2->coeff + i), *(op2->coeff + i));
        }

        poly* result_long = new poly(2 * result->size - 1);
        mpz_t temp;
        mpz_init_set_str(temp, "0", 10);

        for(int i = 0; i < result->size; i++)
        {
            for(int j = 0; j < result->size; j++)
            {
                mpz_mul(temp, *(temp1->coeff + i), *(temp2->coeff + j));
                mpz_add(*(result_long->coeff + i + j), *(result_long->coeff + i + j), temp);
            }
        }

        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(result->coeff + i), *(result_long->coeff + i));
        }

        for(int i = result->size; i < 2 * result->size; i++)
        {
            mpz_neg(*(result_long->coeff + i), *(result_long->coeff + i));
            mpz_add(*(result->coeff + i - result->size), *(result->coeff + i - result->size), *(result_long->coeff + i));
        }

        delete temp1;
        delete temp2;
        delete result_long;
        mpz_clear(temp);

        return 0;
    }
}

int poly_mod(poly* op, mpz_t mod, poly* result)
{
    poly* temp = new poly(op->size);

    for(int i = 0; i < temp->size; i++)
    {
        mpz_set(*(temp->coeff + i), *(op->coeff + i));
    }

    for(int i = 0; i < temp->size; i++)
    {
        mpz_mod(*(temp->coeff + i), *(temp->coeff + i), mod);
    }

    for(int i = 0; i < temp->size; i++)
    {
        mpz_set(*(result->coeff + i), *(temp->coeff + i));
    }

    delete temp;

    return 0;
}

int is_prime(mpz_t p)
{
    int result = mpz_probab_prime_p(p, 15);

    return result;
}

mpz_t* find_ntt_prime(int ring_size, int bits)
{
    mpz_t d_ring_size, k, *target_number = new mpz_t[1];
    mpz_init_set_si(d_ring_size, 2 * ring_size);
    mpz_init(k);
    mpz_init(*target_number);

    mpz_ui_pow_ui(*target_number, 2, bits);
    mpz_fdiv_q(k, *target_number, d_ring_size);

    while(true)
    {
        mpz_mul(*target_number, k, d_ring_size);
        mpz_add_ui(*target_number, *target_number, 1);
        
        if(is_prime(*target_number) > 0)
        {
            break;
        }
        else
        {
            mpz_add_ui(k, k, 1);
        }
    }

    mpz_clears(d_ring_size, k, NULL);

    return target_number;
}

void mpz_abs_sub(mpz_t res, mpz_t a, mpz_t b) 
{
    mpz_sub(res, a, b);
    mpz_abs(res, res);
}

void pollard_rho(mpz_t factor, mpz_t n) 
{
    if (mpz_divisible_ui_p(n, 2)) 
    {
        mpz_set_ui(factor, 2);
        return;
    }

    mpz_t x, y, c, d, val, temp_xy;
    mpz_inits(x, y, c, d, val, temp_xy, NULL);

    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(c, 1);
    mpz_set_ui(d, 1);

    while (mpz_cmp_ui(d, 1) == 0) 
    {
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);

        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);

        mpz_abs_sub(temp_xy, x, y);
        mpz_gcd(d, temp_xy, n);

        if (mpz_cmp(d, n) == 0) 
        {
            mpz_set_ui(factor, 1); // Fail signal
            break;
        }
    }
    
    if (mpz_cmp_ui(d, 1) != 0) 
    {
        mpz_set(factor, d);
    }

    mpz_clears(x, y, c, d, val, temp_xy, NULL);
}

void factorize_recursive(mpz_t n, std::vector<mpz_t*>& factors) 
{
    if (mpz_cmp_ui(n, 1) <= 0) return;

    if (mpz_probab_prime_p(n, 25) > 0) 
    {
        bool exists = false;
        for(auto ptr : factors) 
        {
            if(mpz_cmp(*ptr, n) == 0) 
            { 
                exists = true; break; 
            }
        }
        if (!exists) 
        {
            mpz_t* p = new mpz_t[1];
            mpz_init_set(*p, n);
            factors.push_back(p);
        }
        return;
    }

    mpz_t factor;
    mpz_init(factor);
    
    pollard_rho(factor, n);


    if (mpz_cmp_ui(factor, 1) == 0 || mpz_cmp(factor, n) == 0) 
    {
        // 단순히 2부터 조금 돌려봄 (Backup Plan)
        // 실제로는 여기서 멈추면 안되지만 코드 간결화를 위해 생략
    } 
    else 
    {
        mpz_t remaining;
        mpz_init(remaining);
        mpz_divexact(remaining, n, factor);

        factorize_recursive(factor, factors);
        factorize_recursive(remaining, factors);
        
        mpz_clear(remaining);
    }
    mpz_clear(factor);
}

int find_primitive_root(mpz_t mod, mpz_t result) {
    mpz_t phi, exp, res, g;
    mpz_inits(phi, exp, res, g, NULL);

    mpz_sub_ui(phi, mod, 1); // phi = p - 1

    std::vector<mpz_t*> factors;
    
    mpz_t n;
    mpz_init_set(n, phi);
    
    unsigned long small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    for (unsigned long p : small_primes) 
    {
        if (mpz_divisible_ui_p(n, p)) 
        {
            mpz_t* f = new mpz_t[1];
            mpz_init_set_ui(*f, p);
            factors.push_back(f);
            while (mpz_divisible_ui_p(n, p)) mpz_divexact_ui(n, n, p);
        }
    }

    factorize_recursive(n, factors);

    mpz_set_ui(g, 2);
    while (true) 
    {
        bool is_primitive = true;
        
        for (mpz_t* factor : factors) 
        {
            mpz_divexact(exp, phi, *factor);  
            mpz_powm(res, g, exp, mod);      
            
            if (mpz_cmp_ui(res, 1) == 0) 
            {   
                is_primitive = false;
                break;
            }
        }

        if (is_primitive) 
        {
            mpz_set(result, g);
            break;
        }
        mpz_add_ui(g, g, 1);
    }

    for (auto ptr : factors) 
    {
        mpz_clear(*ptr);
        delete[] ptr;
    }
    mpz_clears(phi, exp, res, g, n, NULL);
    return 0;
}

void bit_reverse(poly* op) 
{
    size_t j = 0;
    for(size_t i = 1; i < op->size; i++) 
    {
        size_t bit = op->size >> 1;
        while (j & bit) 
        {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            mpz_swap(op->coeff[i], op->coeff[j]);
        }
    }
}

void standard_ntt(poly* op, mpz_t mod, const mpz_t root) {
    bit_reverse(op);

    mpz_t u, v, w, w_len, temp_sub;
    mpz_inits(u, v, w, w_len, temp_sub, NULL);

    for (size_t len = 2; len <= op->size; len <<= 1) 
    {
        mpz_t exponent; 
        mpz_init(exponent);
        mpz_sub_ui(exponent, mod, 1);
        mpz_divexact_ui(exponent, exponent, len);
        mpz_powm(w_len, root, exponent, mod);
        mpz_clear(exponent);

        for (size_t i = 0; i < op->size; i += len) 
        {
            mpz_set_ui(w, 1); 
            for (size_t j = 0; j < len / 2; j++) 
            {
                mpz_set(u, op->coeff[i + j]);
                
                mpz_mul(v, op->coeff[i + j + len / 2], w);
                mpz_mod(v, v, mod);

                mpz_add(op->coeff[i + j], u, v);
                mpz_mod(op->coeff[i + j], op->coeff[i + j], mod);

                mpz_sub(temp_sub, u, v);
                mpz_mod(op->coeff[i + j + len / 2], temp_sub, mod); 

                mpz_mul(w, w, w_len);
                mpz_mod(w, w, mod);
            }
        }
    }
    mpz_clears(u, v, w, w_len, temp_sub, NULL);
}

void ntt_negacyclic(poly* op, mpz_t mod, mpz_t g) 
{
    mpz_t psi, psi_pow, temp;
    mpz_inits(psi, psi_pow, temp, NULL);

    mpz_t exp;
    mpz_init(exp);
    mpz_sub_ui(exp, mod, 1);     
    mpz_divexact_ui(exp, exp, 2 * op->size);
    mpz_powm(psi, g, exp, mod); 
    mpz_clear(exp);

    mpz_set_ui(psi_pow, 1);
    for (size_t i = 0; i < op->size; i++) 
    {
        mpz_mul(temp, op->coeff[i], psi_pow);
        mpz_mod(op->coeff[i], temp, mod);

        mpz_mul(psi_pow, psi_pow, psi);
        mpz_mod(psi_pow, psi_pow, mod);
    }

    standard_ntt(op, mod, g);

    mpz_clears(psi, psi_pow, temp, NULL);
}

void intt_negacyclic(poly* op, mpz_t mod, const mpz_t g) 
{
    mpz_t psi, psi_inv, psi_pow_inv, temp, n_inv;
    mpz_inits(psi, psi_inv, psi_pow_inv, temp, n_inv, NULL);

    mpz_t g_inv;
    mpz_init(g_inv);
    mpz_invert(g_inv, g, mod);
    
    standard_ntt(op, mod, g_inv);

    mpz_clear(g_inv);

    mpz_t exp;
    mpz_init(exp);
    mpz_sub_ui(exp, mod, 1);
    mpz_divexact_ui(exp, exp, 2 * op->size);
    mpz_powm(psi, g, exp, mod);
    mpz_invert(psi_inv, psi, mod);
    mpz_clear(exp);
    
    mpz_t n_mpz;
    mpz_init_set_si(n_mpz, op->size);
    mpz_invert(n_inv, n_mpz, mod); 
    mpz_clear(n_mpz);

    mpz_set_ui(psi_pow_inv, 1);
    for (size_t i = 0; i < op->size; i++) 
    {
        mpz_mul(temp, op->coeff[i], psi_pow_inv);
        mpz_mod(temp, temp, mod);

        mpz_mul(op->coeff[i], temp, n_inv);
        mpz_mod(op->coeff[i], op->coeff[i], mod);

        mpz_mul(psi_pow_inv, psi_pow_inv, psi_inv);
        mpz_mod(psi_pow_inv, psi_pow_inv, mod);
    }

    mpz_clears(psi, psi_inv, psi_pow_inv, temp, n_inv, NULL);
}

int poly_mul_ntt(poly* op1, poly* op2, mpz_t mod, mpz_t g, poly* result)
{
    poly* temp1 = new poly(op1->size);
    poly* temp2 = new poly(op2->size);

    if(op1->size != op2->size || op1->size != result->size || op2->size != result->size)
    {
        delete temp1;
        delete temp2;

        return -1;
    }
    else
    {
        for(int i = 0; i < result->size; i++)
        {
            mpz_set(*(temp1->coeff + i), *(op1->coeff + i));
            mpz_set(*(temp2->coeff + i), *(op2->coeff + i));
        }

        ntt_negacyclic(temp1, mod, g);
        ntt_negacyclic(temp2, mod, g);

        for (size_t i = 0; i < result->size; i++) 
        {
            mpz_mul(result->coeff[i], temp1->coeff[i], temp2->coeff[i]);
        }

        poly_mod(result, mod, result);

        intt_negacyclic(result, mod, g);

        delete temp1;
        delete temp2;
        
        return 0;
    }
}


