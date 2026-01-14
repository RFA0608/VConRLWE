#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <gmpxx.h>

class poly
{
    public:
        std::vector<mpz_class> coeff;
        int ring_dim;

    poly(int ring_dim)
    {
        this->coeff.resize(ring_dim);
        this->ring_dim = ring_dim;
    }

    ~poly()
    {
        this->coeff.clear();
    }

    poly* clone()
    {
        poly* clone = new poly(this->ring_dim);
        clone->coeff = this->coeff;

        return clone;
    }
};

class seed_gen
{
    public:
        static gmp_randstate_t& get_rand_state()
        {
            thread_local static struct randwrapper
            {
                gmp_randstate_t state;

                randwrapper()
                {
                    gmp_randinit_default(state);
            
                    std::random_device rd;
                    gmp_randseed_ui(state, rd());
                }

                ~randwrapper()
                {
                    gmp_randclear(state);
                }
            } wrapper;

            return wrapper.state;
        }
};

class poly_handler
{
    public:
        static int plain_2_poly(poly* op, poly* res)
        {
            if(op->ring_dim != res->ring_dim)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = op->coeff[i];
                }

                return 0;
            }
        }

        static int poly_2_plain(poly* op, poly* res)
        {
            if(op->ring_dim != res->ring_dim)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = op->coeff[i];
                }

                return 0;
            }
        }

        static int poly_neg(poly* op, poly* res)
        {
            if(op->ring_dim != res->ring_dim)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = -op->coeff[i];
                }

                return 0;
            }
        }

        static int poly_add(poly* op1, poly* op2, poly* res)
        {
            if((op1->ring_dim != res->ring_dim) || (op2->ring_dim != res->ring_dim))
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = op1->coeff[i] + op2->coeff[i];
                }

                return 0;
            }
        }

        static int poly_scal_mul(mpz_class c, poly* op, poly* res)
        {
            if((op->ring_dim != res->ring_dim))
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = c * op->coeff[i];
                }

                return 0;
            }
        }

        static int poly_mul(poly* op1, poly* op2, poly* res)
        {
            if((op1->ring_dim != res->ring_dim) || (op2->ring_dim != res->ring_dim))
            {
                return -1;
            }
            else
            {
                poly* longres = new poly(op1->ring_dim + op2->ring_dim);
                mpz_class temp;

                for(int i = 0; i < res->ring_dim; i++)
                {
                    for(int j = 0; j < res->ring_dim; j++)
                    {
                        temp = op1->coeff[i] * op2->coeff[j];
                        longres->coeff[i+j] = longres->coeff[i+j] + temp;
                    }
                }

                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = longres->coeff[i];
                }

                for(int i = res->ring_dim; i < 2 * res->ring_dim; i++)
                {
                    res->coeff[i - res->ring_dim] = res->coeff[i - res->ring_dim] - longres->coeff[i];
                }

                delete longres;
                return 0;
            }
        }

        static int poly_mod(poly* op, const mpz_class& mod, poly* res)
        {
            if(op->ring_dim != res->ring_dim)
            {
                return -1;
            }
            else
            {
                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = ((op->coeff[i]) % mod + mod) % mod;
                }

                return 0;
            }
        }

        static int poly_bias_mod(poly* op, const mpz_class& mod, poly* res)
        {
            if(op->ring_dim != res->ring_dim)
            {
                return -1;
            }
            else
            {
                mpz_class half_mod = mod / 2;

                for(int i = 0; i < res->ring_dim; i++)
                {
                    if(op->coeff[i] > half_mod)
                    {
                        res->coeff[i] = ((op->coeff[i] - mod) % mod);
                    }
                    else
                    {
                        res->coeff[i] = (op->coeff[i] % mod);
                    }
                }

                return 0;
            }
        }

        static void bit_reverse(poly* opres)
        {
            size_t j = 0;
            for(size_t i = 1; i < opres->ring_dim; i++)
            {
                size_t bit = opres->ring_dim >> 1;
                while(j & bit)
                {
                    j ^= bit;
                    bit >>= 1;
                }
                j^= bit;
                if(i < j)
                {
                    std::swap(opres->coeff[i], opres->coeff[j]);
                }
            }
        }

        static void standard_ntt_g(poly* op, const mpz_class& q, const mpz_class& g)
        {
            bit_reverse(op);

            mpz_class u, v, w, w_len;

            for(int len = 2; len <= op->ring_dim; len <<= 1)
            {
                mpz_class exp = q - 1;
                exp = exp / len;
                mpz_powm(w_len.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());

                for(int i = 0; i < op->ring_dim; i += len)
                {
                    w = 1;
                    for(int j = 0; j < len / 2; j++)
                    {
                        u = op->coeff[i + j];
                        v = (op->coeff[i + j + len / 2] * w) % q;
                        
                        op->coeff[i + j] = (u + v) % q;
                        op->coeff[i + j + len / 2] = (u - v + q) % q;

                        w = (w * w_len) % q;
                    }
                }
            }
        }

        static void standard_ntt(poly* op, const mpz_class& q, const mpz_class& omega)
        {
            bit_reverse(op);

            mpz_class u, v, w, w_len;

            for(int len = 2; len <= op->ring_dim; len <<= 1)
            {
                mpz_class exp = op->ring_dim;
                exp = exp / len;
                mpz_powm(w_len.get_mpz_t(), omega.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());

                for(int i = 0; i < op->ring_dim; i += len)
                {
                    w = 1;
                    for(int j = 0; j < len / 2; j++)
                    {
                        u = op->coeff[i + j];
                        v = (op->coeff[i + j + len / 2] * w) % q;
                        
                        op->coeff[i + j] = (u + v) % q;
                        op->coeff[i + j + len / 2] = (u - v + q) % q;

                        w = (w * w_len) % q;
                    }
                }
            }
        }

        static void negacyclic_ntt_g(poly* op, const mpz_class& q, const mpz_class& g)
        {
            mpz_class exp = (q - 1) / (2 * op->ring_dim);
            mpz_class psi;
            mpz_powm(psi.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());
            
            mpz_class psi_pow = 1;
            for(int i = 0; i < op->ring_dim; i++)
            {
                op->coeff[i] = (op->coeff[i] * psi_pow) % q;
                psi_pow = (psi_pow * psi) % q;
            }

            standard_ntt(op, q, g);
        }

        static void negacyclic_ntt(poly* op, const mpz_class& q, const mpz_class& psi)
        {
            mpz_class psi_pow = 1;
            for(int i = 0; i < op->ring_dim; i++)
            {
                op->coeff[i] = (op->coeff[i] * psi_pow) % q;
                psi_pow = (psi_pow * psi) % q;
            }

            mpz_class omega;
            mpz_powm_ui(omega.get_mpz_t(), psi.get_mpz_t(), 2, q.get_mpz_t());

            standard_ntt(op, q, omega);
        }

        static void negacyclic_intt_g(poly* op, const mpz_class& q, const mpz_class& g)
        {
            mpz_class g_inv;
            mpz_invert(g_inv.get_mpz_t(), g.get_mpz_t(), q.get_mpz_t());
            
            standard_ntt(op, q, g_inv);

            mpz_class exp = (q - 1) / (2 * op->ring_dim);
            mpz_class psi, psi_inv;
            mpz_powm(psi.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(), q.get_mpz_t());
            mpz_invert(psi_inv.get_mpz_t(), psi.get_mpz_t(), q.get_mpz_t());

            mpz_class n = op->ring_dim;
            mpz_class n_inv;
            mpz_invert(n_inv.get_mpz_t(), n.get_mpz_t(), q.get_mpz_t());

            mpz_class psi_pow_inv = 1;
            for(int i = 0; i < op->ring_dim; i++)
            {
                op->coeff[i] = (n_inv * ((op->coeff[i] * psi_pow_inv) % q)) % q;
                psi_pow_inv = (psi_pow_inv * psi_inv) % q;
            }
        }

        static void negacyclic_intt(poly* op, const mpz_class& q, const mpz_class& psi)
        {
            mpz_class psi_inv, omega_inv;
            mpz_invert(psi_inv.get_mpz_t(), psi.get_mpz_t(), q.get_mpz_t());
            mpz_powm_ui(omega_inv.get_mpz_t(), psi_inv.get_mpz_t(), 2, q.get_mpz_t());

            standard_ntt(op, q, omega_inv);

            mpz_class n = op->ring_dim;
            mpz_class n_inv;
            mpz_invert(n_inv.get_mpz_t(), n.get_mpz_t(), q.get_mpz_t());

            mpz_class psi_pow_inv = 1;
            for(int i = 0; i < op->ring_dim; i++)
            {
                op->coeff[i] = (n_inv * ((op->coeff[i] * psi_pow_inv) % q)) % q;
                psi_pow_inv = (psi_pow_inv * psi_inv) % q;
            }
        }

        static int poly_mul_ntt_g(poly* op1, poly* op2, const mpz_class& mod, const mpz_class& g, poly* res)
        {
            if((op1->ring_dim != res->ring_dim) || (op2->ring_dim != res->ring_dim))
            {
                return -1;
            }
            else
            {
                poly* clone1 = op1->clone();
                poly* clone2 = op2->clone();

                negacyclic_ntt(clone1, mod, g);
                negacyclic_ntt(clone2, mod, g);

                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = (clone1->coeff[i] * clone2->coeff[i]) % mod;
                }

                negacyclic_intt(res, mod, g);


                delete clone1;
                delete clone2;
                return 0;
            }
        }

        static int poly_mul_ntt(poly* op1, poly* op2, const mpz_class& mod, const mpz_class& psi, poly* res)
        {
            if((op1->ring_dim != res->ring_dim) || (op2->ring_dim != res->ring_dim))
            {
                return -1;
            }
            else
            {
                poly* clone1 = op1->clone();
                poly* clone2 = op2->clone();

                negacyclic_ntt(clone1, mod, psi);
                negacyclic_ntt(clone2, mod, psi);

                for(int i = 0; i < res->ring_dim; i++)
                {
                    res->coeff[i] = (clone1->coeff[i] * clone2->coeff[i]) % mod;
                }

                negacyclic_intt(res, mod, psi);

                delete clone1;
                delete clone2;
                return 0;
            }
        }
};

class batch_encoder
{
    public:
        static void encode_g(std::vector<int64_t>& slots, const mpz_class& p, const mpz_class& g, poly* res) 
        {
            for (int i = 0; i < res->ring_dim; i++) 
            {
                if (i < slots.size()) 
                {
                    res->coeff[i] = slots[i];
                } 
                else 
                {
                    res->coeff[i] = 0;
                }
            }

            poly_handler::negacyclic_intt_g(res, p, g);
        }

        static void encode(std::vector<int64_t>& slots, const mpz_class& p, const mpz_class& psi, poly* res) 
        {
            for (int i = 0; i < res->ring_dim; i++) 
            {
                if (i < slots.size()) 
                {
                    res->coeff[i] = slots[i];
                } 
                else 
                {
                    res->coeff[i] = 0;
                }
            }

            poly_handler::negacyclic_intt(res, p, psi);
        }

        static void decode_g(poly* op, const mpz_class& p, const mpz_class& g, std::vector<int64_t>& res) 
        {
            poly* temp = op->clone();

            poly_handler::negacyclic_ntt_g(temp, p, g);

            res.resize(op->ring_dim);

            for (size_t i = 0; i < op->ring_dim; i++) 
            {
                res[i] = temp->coeff[i].get_ui(); 
            }

            delete temp;
        }

        static void decode(poly* op, const mpz_class& p, const mpz_class& psi, std::vector<int64_t>& res) 
        {
            poly* temp = op->clone();

            poly_handler::negacyclic_ntt(temp, p, psi);

            res.resize(op->ring_dim);

            for (size_t i = 0; i < op->ring_dim; i++) 
            {
                res[i] = temp->coeff[i].get_ui(); 
            }

            delete temp;
        }
};

class prime_handler
{
    public:
        static int is_prime(const mpz_class& p)
        {
            int result = mpz_probab_prime_p(p.get_mpz_t(), 25);

            return result;
        }

        static int find_ntt_prime(int ring_dim, int bits, mpz_class& res)
        {
            mpz_class ring_dim_tw = 2 * ring_dim;
            res = mpz_class(1) << bits;
            mpz_class k = res / ring_dim_tw;

            while(true)
            {
                res = k * ring_dim_tw;
                res = res + 1;

                if(is_prime(res) > 0)
                {
                    break;
                }
                else
                {
                    k = k + 1;
                }
            }

            
            return 0;
        }

        static int pollard_rho(const mpz_class& n, mpz_class& res)
        {
            if(n % 2 == 0)
            {
                res = 2;

                return 0;
            }
            else if(is_prime(n) > 0)
            {
                res = n;
                return 0;
            }
            else
            {
                mpz_class x = 2;
                mpz_class y = 2;
                mpz_class c = 1;
                mpz_class d = 1;
                mpz_class diff;

                auto f = [&](const mpz_class& input) -> mpz_class
                {
                    return (input * input + c) % n;
                };

                while(true)
                {
                    x = 2;
                    y = 2;
                    d = 1;

                    while(d == 1)
                    {
                        x = f(x);
                        y = f(f(y));

                        diff = x - y;
                        if(diff < 0)
                        {
                            diff = - diff;
                        }

                        mpz_gcd(d.get_mpz_t(), diff.get_mpz_t(), n.get_mpz_t());
                    }

                    if(d == n)
                    {
                        c = c + 1;
                    }
                    else
                    {
                        break;
                    }
                }
                
                res = d;

                return 0;
            }
        }

        static int factorize_recursive(const mpz_class& n, std::vector<mpz_class>& factors)
        {
            if(n <= 1)
            {
                return -1;
            }
            
            mpz_class in_n = n;

            unsigned long small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
            for(unsigned long p : small_primes)
            {
                while(in_n % p == 0)
                {
                    factors.push_back(p);
                    in_n = in_n / p;
                }

                if(in_n == 1) 
                {
                    return 0;
                }
            }

            if(is_prime(in_n) > 0)
            {
                factors.push_back(in_n);
                return 0;
            }

            mpz_class factor;
            mpz_class remain;

            pollard_rho(in_n, factor);

            if(factor == 1 || factor == in_n)
            {
                factors.push_back(in_n);
                return 0;
            }
            else
            {
                remain = in_n / factor;
                factorize_recursive(factor, factors);
                factorize_recursive(remain, factors);
            }

            return 0;
            
        }

        static int find_primitive_root(const mpz_class& p, mpz_class& res)
        {
            mpz_class phi = p - 1;
            mpz_class exp, com, g = 2;
            std::vector<mpz_class> factors;

            factorize_recursive(phi, factors);

            if (!factors.empty()) 
            {
                std::sort(factors.begin(), factors.end());
                auto last = std::unique(factors.begin(), factors.end());
                factors.erase(last, factors.end());
            }

            while(true)
            {
                bool is_primitive = true;

                for(const auto& factor : factors)
                {
                    exp = phi / factor;
                    mpz_powm(com.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
                    
                    if(com == 1)
                    {
                        is_primitive = false;
                        break;
                    }
                }

                if(is_primitive)
                {
                    res = g;
                    break;
                }
                else
                {
                    g = g + 1;
                }
            }

            return 0;
        }

        static int find_ntt_root(int ring_dim, const mpz_class& p, mpz_class& res)
        {
            mpz_class ring_dim_tw = 2 * ring_dim;
            mpz_class e = (p - 1) / ring_dim_tw;
            mpz_class r = p - 2;
            mpz_class a;
            mpz_class check;

            while(true)
            {
                mpz_urandomm(a.get_mpz_t(), seed_gen::get_rand_state(), r.get_mpz_t());
                a = a + 2;

                mpz_powm(res.get_mpz_t(), a.get_mpz_t(), e.get_mpz_t(), p.get_mpz_t());

                if(res == 1)
                {
                    continue;
                }

                mpz_powm_ui(check.get_mpz_t(), res.get_mpz_t(), ring_dim, p.get_mpz_t());

                if(check == p - 1)
                {
                    break;
                }
            }

            return 0;
        }

        static int find_schnorr_prime(const mpz_class& q, const int& p_bits, mpz_class& res)
        {
            int k_bits = p_bits - mpz_sizeinbase(q.get_mpz_t(), 2) - 1;

            while(true)
            {
                mpz_urandomb(res.get_mpz_t(), seed_gen::get_rand_state(), k_bits);
            
                res = q * res;
                mpz_mul_2exp(res.get_mpz_t(), res.get_mpz_t(), 1);
                res = res + 1;
                
                if(is_prime(res) > 0)
                {
                    break;
                }
            }
            
            return 0;
        }
        
        static int find_schnorr_gen(const mpz_class& q, const mpz_class& p, mpz_class& res)
        {
            mpz_class e = (p - q) / q;
            mpz_class h = 2;

            while(true)
            {
                mpz_powm(res.get_mpz_t(), h.get_mpz_t(), e.get_mpz_t(), p.get_mpz_t());

                if(res != 1)
                {
                    break;
                }
                else
                {
                    h += 1;
                }
            }

            return 0;
        }
};