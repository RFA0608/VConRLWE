#include <random>
#include <ctime>

#include "Struct.h"

class cipher
{
    public:
        std::vector<poly*> ciphertext;
        int level = 1;
        mpz_t plain_mod;
        mpz_t cipher_mod;
        mpz_t g_plain;
        mpz_t g_cipher;
        
        cipher(int size, mpz_t plain_mod, mpz_t cipher_mod)
        {
            mpz_init_set(this->plain_mod, plain_mod);
            mpz_init_set(this->cipher_mod, cipher_mod);
            mpz_init(this->g_plain);
            mpz_init(this->g_cipher);

            find_primitive_root(plain_mod, this->g_plain);
            find_primitive_root(cipher_mod, this->g_cipher);
        }

        ~cipher()
        {
            for(int i = 0; i <= this->level; i++)
            {
                delete this->ciphertext[i];
            }

            mpz_clears(this->plain_mod, this->cipher_mod, this->g_plain, this->g_cipher, NULL);
        }
};

mpz_t* sample_cbd(std::mt19937& gen)
{
    // std::mt19937 gen(std::random_device{}());

    std::uniform_int_distribution<int> dist(0, 1);

    int a = dist(gen) + dist(gen);
    int b = dist(gen) + dist(gen);

    mpz_t* result = new mpz_t[1];
    mpz_init_set_si(*result, a - b);
    
    return result;
}

mpz_t* sample_uni(mpz_t q, gmp_randstate_t state)
{
    // gmp_randstate_t state;
    // gmp_randinit_default(state);
    // gmp_randseed_ui(state, time(NULL));

    mpz_t *result = new mpz_t[1];
    mpz_init(*result);
    mpz_urandomm(*result, state, q);

    return result;
}

poly* secret_key(int size)
{
    poly* key = new poly(size);
    mpz_t *temp = new mpz_t[1];

    std::mt19937 gen(std::random_device{}());
    
    for(int i = 0; i < key->size; i++)
    {
        temp = sample_cbd(gen);
        mpz_set(key->coeff[i], *temp);
        mpz_clear(*temp);
        delete[] temp;
    }
    
    return key;
}

poly* public_key(int size, mpz_t q)
{
    poly* key = new poly(size);
    mpz_t *temp = new mpz_t[1];

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    
    for(int i = 0; i < key->size; i++)
    {
        temp = sample_uni(q, state);
        mpz_set(key->coeff[i], *temp);
        mpz_clear(*temp);
        delete[] temp;
    }
    
    gmp_randclear(state);
    return key;
}

poly* errors(int size, mpz_t p)
{
    poly* err = new poly(size);
    mpz_t *temp = new mpz_t[1];

    std::mt19937 gen(std::random_device{}());
    
    for(int i = 0; i < err->size; i++)
    {
        temp = sample_cbd(gen);
        mpz_set(err->coeff[i], *temp);
        mpz_mul(err->coeff[i], err->coeff[i], p);
        mpz_clear(*temp);
        delete[] temp;
    }
    
    return err;
}

void encrypt(poly* m, poly* sk, cipher* res)
{
    poly* a = public_key(m->size, res->cipher_mod);
    poly* ska = new poly(m->size);

    poly_mul_ntt(sk, a, res->cipher_mod, res->g_cipher, ska);
    
    poly* e = errors(m->size, res->plain_mod);

    poly* b = new poly(m->size);
    poly_add(m, ska, b);
    poly_add(b, e, b);
    poly_mod(b, res->cipher_mod, b);
    
    res->ciphertext.push_back(b);
    res->ciphertext.push_back(a);

    delete ska;
    delete e;
}

int decrypt(cipher* c, poly* sk, poly* res)
{
    switch (c->level)
    {
        case 1:
        {
            poly* ska = new poly(res->size);

            poly_mul_ntt(sk, c->ciphertext[1], c->cipher_mod, c->g_cipher, ska);

            poly_neg(ska, ska);
            poly_add(c->ciphertext[0], ska, res);
            poly_mod(res, c->cipher_mod, res);

            mpz_t half_q;
            mpz_init(half_q);
            mpz_div_ui(half_q, c->cipher_mod, 2);

            for(int i = 0; i < res->size; i++) {
                if(mpz_cmp(res->coeff[i], half_q) > 0) 
                {
                    mpz_sub(res->coeff[i], res->coeff[i], c->cipher_mod);
                }
            }
            mpz_clear(half_q);

            poly_mod(res, c->plain_mod, res);

            delete ska;
            break;
        }
        case 2:
        {

            
            break;
        }
        default:
        {

            return -1;
            break;
        }
    }

    return 0;
}