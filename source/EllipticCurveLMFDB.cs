using System;
using System.Globalization;
using System.Net.Http;
using System.Numerics;
using System.Text;
using System.Text.Json;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveLMFDB
    {
        public EllipticCurveLMFDB(EllipticCurveQ ellipticCurve) 
        {
            GetLmfdbRecord(ellipticCurve);
        }

        #region Properties

        public int Rank
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                return _lmfdbCache.Value.rank;
            }
        }

        public int? AnalyticRank
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                return _lmfdbCache.Value.analyticRank;
            }
        }

        public BigInteger Conductor
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                return _lmfdbCache.Value.conductor;
            }
        }

        public string Label
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                return _lmfdbCache.Value.label;
            }
        }

        public string Url
        {
            get
            {
                if (string.IsNullOrEmpty(Label))
                {
                    return string.Empty;
                }

                return $"https://www.lmfdb.org/EllipticCurve/Q/{Label}/";
            }
        }

        public string TorsionStructure
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                return _lmfdbCache.Value.torsionStructure;
            }
        }

        public EllipticCurveQ GlobalMinimalModel
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                var a = _lmfdbCache.Value.ainvs;
                return new EllipticCurveQ(
                    new BigRational(a[0]),
                    new BigRational(a[1]),
                    new BigRational(a[2]),
                    new BigRational(a[3]),
                    new BigRational(a[4]));
            }
        }

        #region Private methods

        // --- LMFDB cache (ainvs, conductor, rank) ---
        private (BigInteger[] ainvs, 
            BigInteger conductor, 
            int rank, 
            int? analyticRank,
            string label,
            string torsionStructure)? _lmfdbCache;

        // --- Core: fetch from LMFDB by exact j-invariant and match Q-isomorphism ---
        private void GetLmfdbRecord(EllipticCurveQ ellipticCurve)
        {
            if (_lmfdbCache.HasValue) return;

            // 1) Build URL by exact rational j-invariant (e.g. "1556068/81")
            var jinv = ellipticCurve.JInvariant; // BigRational (exact)
            string jnum = jinv.Num.ToString(CultureInfo.InvariantCulture);
            string jden = jinv.Den.ToString(CultureInfo.InvariantCulture);
            string url =
                $"https://www.lmfdb.org/api/ec_curvedata/?jinv=li{jnum},{jden}" +
                "&_format=json" +
                "&_fields=ainvs,conductor,rank,analytic_rank,iso_nlabel,torsion_structure,lmfdb_label";

            // 2) Query LMFDB
            using var httpClient = new HttpClient();
            httpClient.Timeout = new TimeSpan(0, 0, 10);
            httpClient.DefaultRequestHeaders.Add(
                "User-Agent",
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0.0.0 Safari/537.36");

            using var response = httpClient.GetAsync(url).GetAwaiter().GetResult();
            var json = response.Content.ReadAsStringAsync().GetAwaiter().GetResult();

            using var doc = JsonDocument.Parse(json);
            if (!doc.RootElement.TryGetProperty("data", out var data) || data.ValueKind != JsonValueKind.Array)
                throw new InvalidOperationException("LMFDB: unexpected response format");

            // Invariants (as rationals)
            var c4E = ellipticCurve.C4;
            var c6E = ellipticCurve.C6;
            var dE = ellipticCurve.Discriminant;

            // 3) Scan candidates and pick Q-isomorphic one
            for (int i = 0; i < data.GetArrayLength(); i++)
            {
                var row = data[i];
                var ainvs = ParseAinvs(row.GetProperty("ainvs"));
                var (c4C, c6C, dC) = InternalMath.InvariantsIntFromAinvs(ainvs);

                if (InternalMath.IsQIsomorphic(c4E, c6E, dE, c4C, c6C, dC, out _))
                {
                    var conductor = ReadBigInteger(row.GetProperty("conductor"));
                    var rank = row.GetProperty("rank").GetInt32();
                    int? analyticRank = row.TryGetProperty("analytic_rank", out var arEl) ? arEl.GetInt32() : null;
                    var lmfdb_label = row.GetProperty("lmfdb_label").GetString();
                    var torsionStructure = FormatTorsionStructure(row.GetProperty("torsion_structure"));

                    _lmfdbCache = (ainvs, conductor, rank, analyticRank, lmfdb_label, torsionStructure);
                    return;
                }
            }

            throw new InvalidOperationException("LMFDB: no Q-isomorphic curve found for this j-invariant");
        }

        // Accepts: JSON array like [2,4] or a ready string like "Z/2Z x Z/4Z"
        private static string FormatTorsionStructure(JsonElement el)
        {
            if (el.ValueKind == JsonValueKind.Array)
            {
                int n = el.GetArrayLength();
                if (n == 0) return "Z/1Z";
                var sb = new StringBuilder();

                for (int i = 0; i < n; i++)
                {
                    if (i > 0) sb.Append(" x ");
                    int k;
                    if (el[i].ValueKind == JsonValueKind.Number) k = el[i].GetInt32();
                    else if (el[i].ValueKind == JsonValueKind.String) k = int.Parse(el[i].GetString(), CultureInfo.InvariantCulture);
                    else k = 1;
                    sb.Append("Z/").Append(k).Append("Z");
                }

                return sb.ToString();
            }
            if (el.ValueKind == JsonValueKind.String)
            {
                var s = el.GetString() ?? string.Empty;
                return string.IsNullOrWhiteSpace(s) ? "Z/1Z" : s;
            }
            return "Z/1Z";
        }

        // ---- Helpers (parsing + invariants + Q-isomorphism check) ----
        private static BigInteger[] ParseAinvs(JsonElement el)
        {
            // LMFDB returns ainvs either as array or as string "[a1,a2,a3,a4,a6]"
            if (el.ValueKind == JsonValueKind.Array)
            {
                var a = new BigInteger[5];
                for (int i = 0; i < 5; i++)
                    a[i] = ReadBigInteger(el[i]);
                return a;
            }

            if (el.ValueKind == JsonValueKind.String)
            {
                var s = el.GetString() ?? throw new FormatException("LMFDB: ainvs is null string");
                s = s.Trim();
                if (s.StartsWith("[")) s = s.Substring(1);
                if (s.EndsWith("]")) s = s.Substring(0, s.Length - 1);

                var rawParts = s.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (rawParts.Length != 5) throw new FormatException("LMFDB: bad ainvs format");

                var a = new BigInteger[5];
                for (int i = 0; i < 5; i++)
                    a[i] = BigInteger.Parse(rawParts[i].Trim(), CultureInfo.InvariantCulture);
                return a;
            }

            throw new FormatException("LMFDB: unexpected ainvs type");
        }

        // ---- Helpers (parse BigInteger) ----
        private static BigInteger ReadBigInteger(JsonElement el)
        {
            return el.ValueKind switch
            {
                JsonValueKind.Number => BigInteger.Parse(el.GetRawText(), CultureInfo.InvariantCulture),
                JsonValueKind.String => BigInteger.Parse(el.GetString()!, CultureInfo.InvariantCulture),
                _ => throw new FormatException("LMFDB: bad integer field")
            };
        }

        #endregion

        #endregion
    }
}
