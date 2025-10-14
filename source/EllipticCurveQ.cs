using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace EllipticCurves
{
    /// <summary>
    /// Elliptic curve over Q in the general Weierstrass form:
    /// y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    /// All coefficients are exact rationals.
    /// </summary>
    public sealed partial class EllipticCurveQ
    {
        public BigRational A1 { get; }

        public BigRational A2 { get; }

        public BigRational A3 { get; }

        public BigRational A4 { get; }

        public BigRational A6 { get; }

        public EllipticCurveQ(BigRational a1, BigRational a2, BigRational a3, BigRational a4, BigRational a6)
        {
            A1 = a1; A2 = a2; A3 = a3; A4 = a4; A6 = a6;
        }

        #region Invariants (b2, b4, b6, b8), c4, c6, discriminant, j-invariant

        public BigRational B2 => A1 * A1 + BigRational.FromInt(4) * A2;

        public BigRational B4 => BigRational.FromInt(2) * A4 + A1 * A3;

        public BigRational B6 => A3 * A3 + BigRational.FromInt(4) * A6;

        public BigRational B8 => A1 * A1 * A6 + BigRational.FromInt(4) * A2 * A6 - A1 * A3 * A4 + A2 * A3 * A3 - A4 * A4;

        public BigRational C4 => B2 * B2 - BigRational.FromInt(24) * B4;

        public BigRational C6 => -B2 * B2 * B2 + BigRational.FromInt(36) * B2 * B4 - BigRational.FromInt(216) * B6;

        public BigRational Discriminant
        {
            get
            {
                // Δ = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
                var t1 = -(B2 * B2) * B8;
                var t2 = -BigRational.FromInt(8) * B4 * B4 * B4;
                var t3 = -BigRational.FromInt(27) * B6 * B6;
                var t4 = BigRational.FromInt(9) * B2 * B4 * B6;
                return t1 + t2 + t3 + t4;
            }
        }

        public BigRational JInvariant
        {
            get
            {
                var delta = Discriminant;

                if (delta.IsZero) throw new InvalidOperationException("Singular curve: discriminant = 0");

                var c4val = C4;
                return c4val * c4val * c4val / delta; // j = c4^3 / Δ
            }
        }

        public bool IsSingular => Discriminant.IsZero;

        #endregion

        #region Curve membership and group law (general Weierstrass formulas)

        public bool IsOnCurve(EllipticCurvePoint P)
        {
            if (P.IsInfinity) return true;
            var x = P.X; var y = P.Y;
            // y^2 + a1*x*y + a3*y  ?=  x^3 + a2*x^2 + a4*x + a6
            var lhs = y * y + A1 * x * y + A3 * y;
            var rhs = x * x * x + A2 * x * x + A4 * x + A6;
            return lhs == rhs;
        }

        public EllipticCurvePoint Negate(EllipticCurvePoint P)
        {
            if (P.IsInfinity) return P;
            // Negation formula for general Weierstrass: if (x,y) is on E, then
            // -(x,y) = (x, -y - a1*x - a3)
            return new EllipticCurvePoint(P.X, -(P.Y + A1 * P.X + A3));
        }

        public EllipticCurvePoint Add(EllipticCurvePoint P, EllipticCurvePoint Q)
        {
            if (P.IsInfinity) return Q;
            if (Q.IsInfinity) return P;

            // P == -Q ?
            if (P.X == Q.X && P.Y + Q.Y + A1 * Q.X + A3 == BigRational.Zero)
                return EllipticCurvePoint.Infinity;

            BigRational lambda;
            if (P.X != Q.X)
            {
                lambda = (Q.Y - P.Y) / (Q.X - P.X);
            }
            else
            {
                var num = BigRational.FromInt(3) * P.X * P.X
                        + BigRational.FromInt(2) * A2 * P.X
                        + A4 - A1 * P.Y;
                var den = BigRational.FromInt(2) * P.Y + A1 * P.X + A3;
                if (den.IsZero) return EllipticCurvePoint.Infinity;
                lambda = num / den;
            }

            var x3 = lambda * lambda + A1 * lambda - A2 - P.X - Q.X;
            var nu = P.Y - lambda * P.X;
            var y3 = -(lambda + A1) * x3 - A3 - nu;

            return new EllipticCurvePoint(x3, y3);
        }

        public EllipticCurvePoint Double(EllipticCurvePoint P) => Add(P, P);

        public EllipticCurvePoint Multiply(EllipticCurvePoint P, BigInteger n)
        {
            if (n.Sign == 0) return EllipticCurvePoint.Infinity;
            if (n.Sign < 0) return Multiply(Negate(P), BigInteger.Negate(n));
            EllipticCurvePoint acc = EllipticCurvePoint.Infinity;
            EllipticCurvePoint basePt = P;
            var k = n;
            while (k > BigInteger.Zero)
            {
                if (!k.IsEven) acc = Add(acc, basePt);
                basePt = Double(basePt);
                k >>= 1;
            }
            return acc;
        }

        #endregion

        #region Short Weierstrass reduction: y^2 = x^3 + A x + B

        public EllipticCurveQ ShortWeierstrass
        {
            get
            {
                // Short Weierstrass model y^2 = x^3 + A x + B returned as an EllipticCurveQ
                // Step 1: y = y' - (a1*x + a3)/2  →  y'^2 = x^3 + a2' x^2 + a4' x + a6'
                var s1 = A1 / BigRational.FromInt(2);
                var s3 = A3 / BigRational.FromInt(2);
                var a2p = A2 + s1 * s1;                              // a2'
                var a4p = A4 + BigRational.FromInt(2) * s1 * s3;     // a4'
                var a6p = A6 + s3 * s3;                              // a6'

                // Step 2: x = X - a2'/3  →  y'^2 = X^3 + A X + B
                var A = a4p - a2p * a2p / BigRational.FromInt(3);
                var B = a6p - a2p * a4p / BigRational.FromInt(3) + BigRational.FromInt(2) * a2p * a2p * a2p / BigRational.FromInt(27);

                // Short Weierstrass curve has a1=a2=a3=0, a4=A, a6=B
                return new EllipticCurveQ(
                    BigRational.Zero,
                    BigRational.Zero,
                    BigRational.Zero,
                    A,
                    B
                );
            }
        }

        public IEnumerable<EllipticCurvePoint> RationalPoints(int numMax, int denMax)
        {
            yield return EllipticCurvePoint.Infinity;

            var two = BigRational.FromInt(2);

            for (BigInteger n = 1; n <= denMax; n++)
            {
                for (BigInteger m = -numMax; m <= numMax; m++)
                {
                    if (BigInteger.GreatestCommonDivisor(BigInteger.Abs(m), n) != BigInteger.One) continue;

                    var x = new BigRational(m, n);

                    // t = a1*x + a3,  S(x) = RHS + (t^2)/4  (RHS = x^3 + a2*x^2 + a4*x + a6)
                    var t = A1 * x + A3;
                    var rhs = x * x * x + A2 * x * x + A4 * x + A6 + t * t / BigRational.FromInt(4);

                    if (BigRational.IsSquare(rhs, out var yp)) // yp = sqrt(S) in Q
                    {
                        var y1 = yp - t / two;
                        var y2 = (-yp) - t / two;

                        var P1 = new EllipticCurvePoint(x, y1);
                        if (IsOnCurve(P1)) yield return P1;

                        if (y2 != y1)
                        {
                            var P2 = new EllipticCurvePoint(x, y2);
                            if (IsOnCurve(P2)) yield return P2;
                        }
                    }
                }
            }
        }

        public IEnumerable<EllipticCurvePoint> IntegralPoints(int xmax)
        {
            return RationalPoints(xmax, 1);
        }

        #endregion

        #region Torsion structure

        private HashSet<EllipticCurvePoint> torsionPoints;

        public IEnumerable<EllipticCurvePoint> TorsionPoints
        {
            get
            {
                if (torsionPoints == null)
                {
                    // ---- 0) Short → integral short model ----
                    var Es = ShortWeierstrass; // y^2 = x^3 + A x + B
                    var A = Es.A4; var B = Es.A6;

                    var d = InternalMath.Lcm(A.Den, B.Den);
                    var Aint = A * BigRational.FromFraction(BigInteger.Pow(d, 4), BigInteger.One);
                    var Bint = B * BigRational.FromFraction(BigInteger.Pow(d, 6), BigInteger.One);
                    if (Aint.Den != BigInteger.One || Bint.Den != BigInteger.One)
                        throw new InvalidOperationException("Failed to obtain integral short model.");

                    var Ashort = Aint.Num;                  // integer A'
                    var Bshort = Bint.Num;                  // integer B'
                    var Eint = new EllipticCurveQ(
                        BigRational.Zero, BigRational.Zero, BigRational.Zero,
                        new BigRational(Ashort), new BigRational(Bshort));

                    // Δ' = -16(4A'^3 + 27B'^2)
                    BigInteger Delta = -16 * (4 * BigInteger.Pow(Ashort, 3) + 27 * BigInteger.Pow(Bshort, 2));
                    if (Delta.IsZero) throw new InvalidOperationException("Singular curve.");

                    // ---- 1) Candidate orders via reductions mod small good primes ----
                    int[] primes = [5, 7, 11, 13, 17, 19, 23, 29];
                    BigInteger gcdOrders = BigInteger.Zero;
                    for (int i = 0; i < primes.Length; i++)
                    {
                        int p = primes[i];
                        if (Delta % p == 0) continue; // good reduction only
                        var np = InternalMath.CountPointsFpShort(Ashort, Bshort, p); // #E(F_p)
                        gcdOrders = gcdOrders.IsZero ? new BigInteger(np) : BigInteger.GreatestCommonDivisor(gcdOrders, np);
                        if (gcdOrders.IsOne) break;
                    }
                    int[] mazur = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
                    var possibleOrders = new List<int>();
                    for (int i = 0; i < mazur.Length; i++)
                    {
                        int n = mazur[i];
                        if (gcdOrders.IsZero || (gcdOrders % n) == 0) possibleOrders.Add(n);
                    }

                    // ---- 2) Lutz–Nagell search on integral short model ----
                    var T = new HashSet<EllipticCurvePoint> { EllipticCurvePoint.Infinity };

                    // 2a) 2-torsion: X | B'  and X^3 + A'X + B' = 0
                    foreach (var x in InternalMath.EnumerateDivisorsAbs(Bshort))
                    {
                        if (InternalMath.EvalCubic(Ashort, Bshort, x) == 0)
                            T.Add(new EllipticCurvePoint(new BigRational(x), BigRational.Zero));
                        if (x != 0 && InternalMath.EvalCubic(Ashort, Bshort, -x) == 0)
                            T.Add(new EllipticCurvePoint(new BigRational(-x), BigRational.Zero));
                    }
                    // Edge: B'=0 ⇒ X = ±sqrt(-A') (if square)
                    if (Bshort.IsZero)
                    {
                        var negA = BigInteger.Negate(Ashort);
                        if (negA.Sign >= 0)
                        {
                            var s = InternalMath.IntegerSqrt(negA);
                            if (s * s == negA && s != 0)
                            {
                                if (InternalMath.EvalCubic(Ashort, Bshort, s) == 0)
                                    T.Add(new EllipticCurvePoint(new BigRational(s), BigRational.Zero));
                                if (InternalMath.EvalCubic(Ashort, Bshort, -s) == 0)
                                    T.Add(new EllipticCurvePoint(new BigRational(-s), BigRational.Zero));
                            }
                        }
                    }

                    // 2b) odd torsion: y != 0 and y^2 | |Δ'|
                    var factDelta = InternalMath.FactorAbs(Delta);
                    foreach (var y2 in InternalMath.EnumerateSquareDivisors(factDelta)) // y2 ≥ 1
                    {
                        if (y2.IsZero) continue;
                        var y = InternalMath.IntegerSqrt(y2); // exact sqrt

                        var C = Bshort - y2;     // X^3 + A'X + (B' - y^2) = 0
                        if (C.IsZero)
                        {
                            var P1 = new EllipticCurvePoint(new BigRational(0), new BigRational(y));
                            var P2 = new EllipticCurvePoint(new BigRational(0), new BigRational(-y));
                            if (Eint.IsOnCurve(P1) && InternalMath.IsTorsionWithCandidates(Eint, P1, possibleOrders)) T.Add(P1);
                            if (Eint.IsOnCurve(P2) && InternalMath.IsTorsionWithCandidates(Eint, P2, possibleOrders)) T.Add(P2);

                            var negA = BigInteger.Negate(Ashort);
                            if (negA.Sign >= 0)
                            {
                                var s = InternalMath.IntegerSqrt(negA);
                                if (s * s == negA && s != 0)
                                {
                                    var xs = new[] { s, BigInteger.Negate(s) };
                                    for (int t = 0; t < xs.Length; t++)
                                    {
                                        var x = xs[t];
                                        if (InternalMath.EvalCubic(Ashort, Bshort, x) == y2)
                                        {
                                            var Q1 = new EllipticCurvePoint(new BigRational(x), new BigRational(y));
                                            var Q2 = new EllipticCurvePoint(new BigRational(x), new BigRational(-y));
                                            if (Eint.IsOnCurve(Q1) && InternalMath.IsTorsionWithCandidates(Eint, Q1, possibleOrders)) T.Add(Q1);
                                            if (Eint.IsOnCurve(Q2) && InternalMath.IsTorsionWithCandidates(Eint, Q2, possibleOrders)) T.Add(Q2);
                                        }
                                    }
                                }
                            }
                            continue;
                        }

                        var absC = C >= 0 ? C : -C;
                        foreach (var x in InternalMath.EnumerateDivisorsAbs(absC))
                        {
                            if (InternalMath.EvalCubic(Ashort, Bshort, x) == y2)
                            {
                                var R1 = new EllipticCurvePoint(new BigRational(x), new BigRational(y));
                                var R2 = new EllipticCurvePoint(new BigRational(x), new BigRational(-y));
                                if (Eint.IsOnCurve(R1) && InternalMath.IsTorsionWithCandidates(Eint, R1, possibleOrders)) T.Add(R1);
                                if (Eint.IsOnCurve(R2) && InternalMath.IsTorsionWithCandidates(Eint, R2, possibleOrders)) T.Add(R2);
                            }
                            if (x != 0)
                            {
                                var nx = -x;
                                if (InternalMath.EvalCubic(Ashort, Bshort, nx) == y2)
                                {
                                    var R1 = new EllipticCurvePoint(new BigRational(nx), new BigRational(y));
                                    var R2 = new EllipticCurvePoint(new BigRational(nx), new BigRational(-y));
                                    if (Eint.IsOnCurve(R1) && InternalMath.IsTorsionWithCandidates(Eint, R1, possibleOrders)) T.Add(R1);
                                    if (Eint.IsOnCurve(R2) && InternalMath.IsTorsionWithCandidates(Eint, R2, possibleOrders)) T.Add(R2);
                                }
                            }
                        }
                    }

                    // ---- 3) Map SHORT-INTEGRAL torsion back to ORIGINAL coordinates ----
                    var mapped = new HashSet<EllipticCurvePoint>
                    {
                        EllipticCurvePoint.Infinity
                    };

                    var d2 = BigInteger.Pow(d, 2);
                    var d3 = BigInteger.Pow(d, 3);

                    // a2' = A2 + (A1/2)^2
                    var s1 = A1 / BigRational.FromInt(2);       // ORIGINAL A1
                    var a2p = A2 + s1 * s1;                     // ORIGINAL A2
                    var shift = a2p / BigRational.FromInt(3);   // a2'/3

                    foreach (var P in T)
                    {
                        if (P.IsInfinity) continue;

                        // short rational coords
                        var X = P.X / new BigRational(d2);
                        var Y = P.Y / new BigRational(d3);

                        // original coords: x = X - a2'/3,  y = Y - (A1*x + A3)/2
                        var x = X - shift;
                        var y = Y - (A1 * x + A3) / BigRational.FromInt(2);

                        var Porig = new EllipticCurvePoint(x, y);
                        if (IsOnCurve(Porig)) mapped.Add(Porig);
                    }

                    torsionPoints = mapped;
                }

                return torsionPoints;
            }
        }

        public string TorsionStructure
        {
            get
            {
                // Deduce group structure from torsion set
                int size = TorsionPoints.Count(); // includes Infinity; equals |E(Q)_tors|
                if (size == 1) return "Z/1Z";

                int twoTors = 0;
                foreach (var P in TorsionPoints)
                {
                    if (!P.IsInfinity && P.Y.IsZero) twoTors++;
                }

                if (twoTors == 3)
                {
                    // E(Q)_tors ≅ Z/2Z × Z/2mZ with m ∈ {1,2,3,4} and |E(Q)_tors| = 4m
                    int m = size / 4;
                    return $"Z/2Z x Z/{2 * m}Z";
                }
                else
                {
                    // cyclic case
                    return $"Z/{size}Z";
                }
            }
        }

        #endregion

        #region Overrides

        public override string ToString()
        {
            var sb = new StringBuilder();

            // ------------ EC ------------
            // Left: y^2 + a1*x*y + a3*y
            sb.Append("y^2");
            AppendTerm(sb, A1, "x*y");
            AppendTerm(sb, A3, "y");

            // Separator
            sb.Append(" = ");

            // Right: x^3 + a2*x^2 + a4*x + a6
            sb.Append("x^3");
            AppendTerm(sb, A2, "x^2");
            AppendTerm(sb, A4, "x");
            AppendTerm(sb, A6, string.Empty);

            return sb.ToString();
        }

        private static void AppendTerm(StringBuilder sb, BigRational coeff, string monomial)
        {
            if (coeff.IsZero) return;

            bool positive = coeff.Sign > 0;
            var abs = positive ? coeff : coeff.Negate();

            sb.Append(positive ? " + " : " - ");

            var showCoeff = string.IsNullOrEmpty(monomial) || abs != BigRational.One;

            if (showCoeff)
            {
                sb.Append(abs.ToString());
                if (!string.IsNullOrEmpty(monomial)) sb.Append('*');
            }

            if (!string.IsNullOrEmpty(monomial))
            {
                sb.Append(monomial);
            }
        }

        #endregion

    }
}
