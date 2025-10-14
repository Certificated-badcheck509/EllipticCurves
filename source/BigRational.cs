using System;
using System.Globalization;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>
    /// Minimal, GC-friendly Big Rational type for exact arithmetic over Q.
    /// Denominator is always > 0 and gcd(|num|, den) = 1.
    /// </summary>
    public readonly struct BigRational : IEquatable<BigRational>, IComparable<BigRational>
    {
        public BigInteger Num { get; }

        public BigInteger Den { get; }

        public static readonly BigRational Zero = new BigRational(BigInteger.Zero, BigInteger.One);

        public static readonly BigRational One = new BigRational(BigInteger.One, BigInteger.One);

        public static readonly BigRational Two = new BigRational(2);

        public BigRational(long n) : this(new BigInteger(n), BigInteger.One) { }

        public BigRational(BigInteger n) : this(n, BigInteger.One) { }

        public BigRational(BigInteger num, BigInteger den)
        {
            if (den.IsZero) throw new DivideByZeroException("Rational with zero denominator");
            if (num.IsZero)
            {
                Num = BigInteger.Zero; Den = BigInteger.One; return;
            }
            if (den.Sign < 0) { num = BigInteger.Negate(num); den = BigInteger.Negate(den); }
            var g = BigInteger.GreatestCommonDivisor(BigInteger.Abs(num), den);
            Num = num / g; Den = den / g;
        }

        public static BigRational FromInt(long n) => new BigRational(n);

        public static BigRational FromFraction(BigInteger num, BigInteger den) => new BigRational(num, den);

        public bool IsZero => Num.IsZero;

        public int Sign => Num.Sign;

        public BigRational Abs() => new BigRational(BigInteger.Abs(Num), Den);

        public BigRational Negate() => new BigRational(BigInteger.Negate(Num), Den);

        public BigRational Reciprocal()
        {
            if (IsZero) throw new DivideByZeroException();
            return new BigRational(Den, Num);
        }

        public static BigRational operator +(BigRational a, BigRational b) => new BigRational(a.Num * b.Den + b.Num * a.Den, a.Den * b.Den);

        public static BigRational operator -(BigRational a, BigRational b) => new BigRational(a.Num * b.Den - b.Num * a.Den, a.Den * b.Den);

        public static BigRational operator *(BigRational a, BigRational b) => new BigRational(a.Num * b.Num, a.Den * b.Den);

        public static BigRational operator /(BigRational a, BigRational b)
        {
            if (b.IsZero) throw new DivideByZeroException();
            return new BigRational(a.Num * b.Den, a.Den * b.Num);
        }

        public static BigRational Pow(BigRational r, int k)
        {
            if (k == 0) return BigRational.One;
            var num = BigInteger.Pow(r.Num, k);
            var den = BigInteger.Pow(r.Den, k);
            return new BigRational(num, den);
        }

        public static bool IsSquare(BigRational r, out BigRational s)
        {
            s = default;
            if (r.Sign < 0) return false;
            var an = BigInteger.Abs(r.Num);
            var ad = r.Den;

            if (!TryIntegerSquareRoot(an, out var sn)) return false;
            if (!TryIntegerSquareRoot(ad, out var sd)) return false;

            s = new BigRational(sn, sd);
            return (s * s) == r; // paranoid check
        }

        private static bool TryIntegerSquareRoot(BigInteger n, out BigInteger rt)
        {
            rt = IntegerSquareRoot(n);
            return rt * rt == n;
        }

        private static BigInteger IntegerSquareRoot(BigInteger n)
        {
            if (n <= 1) return n;
            BigInteger x0 = n, x1 = (n >> 1) + 1;
            while (x1 < x0) { x0 = x1; x1 = (x1 + n / x1) >> 1; }
            return x0;
        }

        public static BigRational operator -(BigRational a) => a.Negate();

        public static bool operator ==(BigRational a, BigRational b) => a.Num == b.Num && a.Den == b.Den;

        public static bool operator !=(BigRational a, BigRational b) => !(a == b);

        public static bool operator <(BigRational a, BigRational b) => a.CompareTo(b) < 0;

        public static bool operator >(BigRational a, BigRational b) => a.CompareTo(b) > 0;

        public static bool operator <=(BigRational a, BigRational b) => a.CompareTo(b) <= 0;

        public static bool operator >=(BigRational a, BigRational b) => a.CompareTo(b) >= 0;

        public int CompareTo(BigRational other)
        {
            var lhs = Num * other.Den;
            var rhs = other.Num * Den;
            return lhs.CompareTo(rhs);
        }

        public bool Equals(BigRational other) => this == other;

        public override bool Equals(object obj) => obj is BigRational r && Equals(r);

        public override int GetHashCode()
        {
            unchecked
            {
                int h = 17;
                h = h * 31 + Num.GetHashCode();
                h = h * 31 + Den.GetHashCode();
                return h;
            }
        }

        public override string ToString()
        {
            if (Den.IsOne) return Num.ToString(CultureInfo.InvariantCulture);
            return $"{Num.ToString(CultureInfo.InvariantCulture)}/{Den.ToString(CultureInfo.InvariantCulture)}";
        }
    }
}
