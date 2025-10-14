using System;

namespace EllipticCurves
{
    /// <summary>
    /// Point on an elliptic curve over Q in affine coords.
    /// </summary>
    public readonly struct EllipticCurvePoint : IEquatable<EllipticCurvePoint>
    {
        public readonly BigRational X;

        public readonly BigRational Y;

        public readonly bool IsInfinity;

        public EllipticCurvePoint(BigRational x, BigRational y) { X = x; Y = y; IsInfinity = false; }

        private EllipticCurvePoint(bool inf) { X = BigRational.Zero; Y = BigRational.Zero; IsInfinity = inf; }

        public static EllipticCurvePoint Infinity => new EllipticCurvePoint(true);

        public bool Equals(EllipticCurvePoint other)
        {
            if (IsInfinity || other.IsInfinity) return IsInfinity == other.IsInfinity;
            return X == other.X && Y == other.Y;
        }

        public override bool Equals(object obj) => obj is EllipticCurvePoint p && Equals(p);

        public override int GetHashCode()
        {
            unchecked
            {
                int h = 17;
                h = h * 31 + X.GetHashCode();
                h = h * 31 + Y.GetHashCode();
                h = h * 31 + (IsInfinity ? 1 : 0);
                return h;
            }
        }

        public override string ToString() => IsInfinity ? "O" : $"({X}, {Y})";
    }
}
