using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    internal static class InternalMath
    {
        // Compute integer c4, c6, Δ for given integral ainvs
        public static (BigInteger c4, BigInteger c6, BigInteger Delta) InvariantsIntFromAinvs(BigInteger[] a)
        {
            var a1 = a[0]; var a2 = a[1]; var a3 = a[2]; var a4 = a[3]; var a6 = a[4];
            BigInteger b2 = a1 * a1 + 4 * a2;
            BigInteger b4 = 2 * a4 + a1 * a3;
            BigInteger b6 = a3 * a3 + 4 * a6;
            BigInteger b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
            BigInteger Delta = -(b2 * b2) * b8 - 8 * b4 * b4 * b4 - 27 * b6 * b6 + 9 * b2 * b4 * b6;
            BigInteger c4 = b2 * b2 - 24 * b4;
            BigInteger c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
            return (c4, c6, Delta);
        }

        // Check c4_E = u^4 * c4_C, c6_E = u^6 * c6_C, Δ_E = u^12 * Δ_C for some u ∈ Q
        public static bool IsQIsomorphic(BigRational c4E, BigRational c6E, BigRational dE, BigInteger c4C, BigInteger c6C, BigInteger dC, out BigRational u)
        {
            u = default;
            var rc4C = new BigRational(c4C);
            var rc6C = new BigRational(c6C);
            var rdC = new BigRational(dC);

            bool c4zero = c4E.IsZero || rc4C.IsZero;
            bool c6zero = c6E.IsZero || rc6C.IsZero;

            if (!c4zero)
            {
                var ratio4 = c4E / rc4C;              // should be u^4
                if (!TryRationalKthRoot(ratio4, 4, out var u4)) return false;
                u = u4;

                if (!c6zero && BigRational.Pow(u, 6) != (c6E / rc6C)) return false;
                if (BigRational.Pow(u, 12) != (dE / rdC)) return false;
                return true;
            }
            else if (!c6zero)
            {
                var ratio6 = c6E / rc6C;              // should be u^6
                if (!TryRationalKthRoot(ratio6, 6, out var u6)) return false;
                u = u6;

                if (BigRational.Pow(u, 12) != (dE / rdC)) return false;
                return true;
            }
            else
            {
                // Both zero cannot happen for a nonsingular curve
                return false;
            }
        }

        // Try to find rational k-th root (k = 4 or 6) of r ∈ Q
        public static bool TryRationalKthRoot(BigRational r, int k, out BigRational root)
        {
            root = default;
            if (r.Sign < 0 && (k % 2 == 0)) return false;

            var an = BigInteger.Abs(r.Num);
            var ad = r.Den;

            if (!TryIntegerKthRoot(an, k, out var rn)) return false;
            if (!TryIntegerKthRoot(ad, k, out var rd)) return false;

            if (r.Sign < 0) rn = BigInteger.Negate(rn);
            root = new BigRational(rn, rd);
            return BigRational.Pow(root, k) == r;
        }

        // Integer k-th root with exactness check
        public static bool TryIntegerKthRoot(BigInteger n, int k, out BigInteger rt)
        {
            rt = IntegerKthRoot(n, k);
            return BigInteger.Pow(rt, k) == n;
        }

        // Integer k-th root
        public static BigInteger IntegerKthRoot(BigInteger n, int k)
        {
            if (n.IsZero) return BigInteger.Zero;
            if (n.IsOne) return BigInteger.One;

            // Rough initial guess: 2^(bitlen/k)
            int bits = n.ToByteArray().Length * 8;
            BigInteger x = BigInteger.One << Math.Max(1, bits / k);

            while (true)
            {
                // Newton iteration: x_{t+1} = ((k-1)*x + n/x^{k-1}) / k
                var xk_1 = BigInteger.Pow(x, k - 1);
                if (xk_1.IsZero) break;
                var next = ((k - 1) * x + n / xk_1) / k;
                if (next >= x) return x; // converged
                x = next;
            }
            return x;
        }


        // Check torsion using candidate orders first; if that fails, fall back to 1..12
        public static bool IsTorsionWithCandidates(EllipticCurveQ E, EllipticCurvePoint P, List<int> orders)
        {
            // Fast path: only divisors of gcd(#E(F_p))
            for (int i = 0; i < orders.Count; i++)
            {
                int n = orders[i];
                if (n <= 1) continue;
                if (OrderDivides(E, P, n)) return true;
            }
            // Fallback: exact check up to Mazur bound (cheap; very few points to test)
            var Q = P;
            for (int k = 1; k <= 12; k++)
            {
                if (Q.IsInfinity) return true; // order divides k
                Q = E.Add(Q, P);
            }
            return false;
        }

        public static bool OrderDivides(EllipticCurveQ E, EllipticCurvePoint P, int n)
        {
            // fast power-like: compute Q = n*P via double-and-add
            var Q = E.Multiply(P, new BigInteger(n));
            return Q.IsInfinity;
        }

        public static BigInteger Lcm(BigInteger a, BigInteger b)
        {
            return a / BigInteger.GreatestCommonDivisor(a, b) * b;
        }

        public static BigInteger EvalCubic(BigInteger A, BigInteger B, BigInteger x)
        {
            return x * x * x + A * x + B;
        }

        public static BigInteger IntegerSqrt(BigInteger n)
        {
            if (n <= 1) return n;
            BigInteger x0 = n, x1 = (n >> 1) + 1;
            while (x1 < x0) { x0 = x1; x1 = (x1 + n / x1) >> 1; }
            return x0;
        }

        // Count points on short Weierstrass mod p (p small prime)
        public static int CountPointsFpShort(BigInteger A, BigInteger B, int p)
        {
            int cnt = 1; // point at infinity
            for (int x = 0; x < p; x++)
            {
                BigInteger rhs = (BigInteger)x * x % p;
                rhs = (rhs * x + A) % p;         // x^3 + A x
                rhs = (rhs + B) % p;             // + B
                if (rhs < 0) rhs += p;

                int chi = Legendre((int)rhs, p); // 0, 1, -1
                if (chi == 0) cnt += 1;          // one solution y=0
                else if (chi == 1) cnt += 2;     // two square roots mod p
            }
            return cnt;
        }

        public static int Legendre(int a, int p)
        {
            if (a % p == 0) return 0;
            // Euler's criterion: a^((p-1)/2) mod p
            int e = (p - 1) / 2;
            int r = ModPow(a, e, p);
            if (r == 1) return 1;
            if (r == p - 1) return -1;
            return 0;
        }

        public static int ModPow(int a, int e, int m)
        {
            long res = 1, b = ((a % m) + m) % m;
            while (e > 0)
            {
                if ((e & 1) == 1) res = (res * b) % m;
                b = b * b % m;
                e >>= 1;
            }
            return (int)res;
        }

        // Factor |n| via Pollard–Rho (enough for our sizes); return p->e
        public static Dictionary<BigInteger, int> FactorAbs(BigInteger n)
        {
            var res = new Dictionary<BigInteger, int>();
            if (n < 0) n = BigInteger.Abs(n);
            if (n <= 1) return res;

            // remove 2s
            int e2 = 0;
            while ((n & 1) == 0) { n >>= 1; e2++; }
            if (e2 > 0) res[new BigInteger(2)] = e2;

            if (n > 1) FactorRec(n, res);
            return res;
        }

        public static void FactorRec(BigInteger n, Dictionary<BigInteger, int> res)
        {
            if (n == 1) return;
            if (IsProbablePrime(n))
            {
                int e;
                if (res.TryGetValue(n, out e)) res[n] = e + 1; else res[n] = 1;
                return;
            }
            var d = PollardRho(n);
            FactorRec(d, res);
            FactorRec(n / d, res);
        }

        public static bool IsProbablePrime(BigInteger n)
        {
            if (n < 2) return false;
            int[] small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
            for (int i = 0; i < small.Length; i++)
            {
                int p = small[i];
                if (n == p) return true;
                if (n % p == 0) return n == p;
            }
            // Miller–Rabin bases good for 64-bit; for larger n ok in practice
            int[] bases = [2, 3, 5, 7, 11, 13, 17];
            BigInteger d = n - 1; int s = 0;
            while ((d & 1) == 0) { d >>= 1; s++; }
            for (int i = 0; i < bases.Length; i++)
            {
                int a = bases[i];
                if (n == a) return true;
                if (MillerRabinCheck(n, d, s, a) == false) return false;
            }
            return true;
        }

        public static bool MillerRabinCheck(BigInteger n, BigInteger d, int s, int a)
        {
            BigInteger x = BigInteger.ModPow(a, d, n);
            if (x == 1 || x == n - 1) return true;
            for (int r = 1; r < s; r++)
            {
                x = x * x % n;
                if (x == n - 1) return true;
            }
            return false;
        }

        public static BigInteger PollardRho(BigInteger n)
        {
            if ((n & 1) == 0) return 2;
            var rnd = new Random(1234567);
            while (true)
            {
                BigInteger c = RandomBelow(n, rnd);
                BigInteger x = RandomBelow(n, rnd);
                BigInteger y = x;
                BigInteger d = 1;
                while (d == 1)
                {
                    x = (x * x + c) % n;
                    y = (y * y + c) % n;
                    y = (y * y + c) % n;
                    d = BigInteger.GreatestCommonDivisor(BigInteger.Abs(x - y), n);
                }
                if (d != n) return d;
            }
        }

        public static BigInteger RandomBelow(BigInteger n, Random rng)
        {
            var bytes = n.ToByteArray();
            BigInteger r;
            do
            {
                rng.NextBytes(bytes);
                bytes[bytes.Length - 1] &= 0x7F;
                r = new BigInteger(bytes);
            } while (r <= 0 || r >= n);
            return r;
        }

        // Enumerate all square divisors y2 of |Δ| using its factorization
        public static IEnumerable<BigInteger> EnumerateSquareDivisors(Dictionary<BigInteger, int> fact)
        {
            // y2 = Π p^{2f}, 0<=f<=floor(e/2)
            var primes = new List<BigInteger>(fact.Keys);
            var exps = new List<int>(primes.Count);
            for (int i = 0; i < primes.Count; i++) exps.Add(fact[primes[i]] / 2);

            var stack = new BigInteger[primes.Count + 1];
            stack[0] = BigInteger.One;

            // DFS explicit to avoid local funcs (older compilers)
            var accs = new List<BigInteger>
            {
                BigInteger.One
            };
            for (int i = 0; i < primes.Count; i++)
            {
                var next = new List<BigInteger>();
                BigInteger p = primes[i];
                int e = exps[i];
                for (int j = 0; j < accs.Count; j++)
                {
                    BigInteger baseVal = accs[j];
                    BigInteger pow = BigInteger.One;
                    for (int t = 0; t <= e; t++)
                    {
                        next.Add(baseVal * pow);
                        pow *= p * p;
                    }
                }
                accs = next;
            }
            for (int i = 0; i < accs.Count; i++) yield return accs[i];
        }

        // Enumerate all positive divisors of |n|
        public static IEnumerable<BigInteger> EnumerateDivisorsAbs(BigInteger n)
        {
            n = n >= 0 ? n : -n;
            if (n.IsZero) { yield return BigInteger.Zero; yield break; }
            var fact = FactorAbs(n);
            var primes = new List<BigInteger>(fact.Keys);
            var exps = new List<int>(primes.Count);
            for (int i = 0; i < primes.Count; i++) exps.Add(fact[primes[i]]);

            var accs = new List<BigInteger>
            {
                BigInteger.One
            };
            for (int i = 0; i < primes.Count; i++)
            {
                var next = new List<BigInteger>();
                BigInteger p = primes[i];
                int e = exps[i];
                for (int j = 0; j < accs.Count; j++)
                {
                    BigInteger baseVal = accs[j];
                    BigInteger pow = BigInteger.One;
                    for (int t = 0; t <= e; t++)
                    {
                        next.Add(baseVal * pow);
                        pow *= p;
                    }
                }
                accs = next;
            }
            for (int i = 0; i < accs.Count; i++) yield return accs[i];
        }
    }
}
