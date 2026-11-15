#ifndef NIS_TESTS_EXACT_BATEMAN_HPP
#define NIS_TESTS_EXACT_BATEMAN_HPP

// The exact solution of a linear decay chain (the Bateman equations), for the solver tests.
//
// The Bateman sum has terms of both signs, so for a small lambda*t it loses many digits. The sum is done in
// double-double arithmetic (a value is hi + lo, two doubles, ~32 significant digits). This uses only double
// and fma, so the result and its error are the same on all platforms. Do not use long double for this: it has
// 64 bits of mantissa on x86-64 Linux, but only 53 (the same as double) on arm64 (macOS on Apple Silicon) and
// on Windows, and then the sum is not precise enough to check the solver to 1e-10.

#include <cmath>
#include <vector>

namespace nis::tests {

    struct DD {
        double hi = 0.0;
        double lo = 0.0;
    };

    namespace dd {

        // a + b = s + e exactly
        inline DD two_sum(double a, double b){
            const double s = a + b;
            const double bb = s - a;
            return {s, (a - (s - bb)) + (b - bb)};
        }

        // a + b = s + e exactly, if |a| >= |b|
        inline DD quick_two_sum(double a, double b){
            const double s = a + b;
            return {s, b - (s - a)};
        }

        // a * b = p + e exactly (with a fused multiply-add)
        inline DD two_prod(double a, double b){
            const double p = a*b;
            return {p, std::fma(a, b, -p)};
        }

        inline DD add(DD a, DD b){
            DD s = two_sum(a.hi, b.hi);
            const DD t = two_sum(a.lo, b.lo);
            s.lo += t.hi;
            s = quick_two_sum(s.hi, s.lo);
            s.lo += t.lo;
            return quick_two_sum(s.hi, s.lo);
        }

        inline DD neg(DD a){
            return {-a.hi, -a.lo};
        }

        inline DD sub(DD a, DD b){
            return add(a, neg(b));
        }

        inline DD mul(DD a, DD b){
            DD p = two_prod(a.hi, b.hi);
            p.lo += a.hi*b.lo + a.lo*b.hi;
            return quick_two_sum(p.hi, p.lo);
        }

        inline DD div(DD a, DD b){
            const double q1 = a.hi/b.hi;
            DD r = sub(a, mul({q1, 0.0}, b));
            const double q2 = r.hi/b.hi;
            r = sub(r, mul({q2, 0.0}, b));
            const double q3 = r.hi/b.hi;
            return add(quick_two_sum(q1, q2), {q3, 0.0});
        }

        inline DD ldexp(DD a, int k){
            return {std::ldexp(a.hi, k), std::ldexp(a.lo, k)};
        }

        // exp(x) for x <= 0, to ~1e-30 relative (less for a result near the smallest double)
        inline DD exp(DD x){
            if(x.hi == 0.0){
                return {1.0, 0.0};
            }
            if(x.hi < -746.0){
                return {0.0, 0.0};
            }
            // x = k*ln(2) + r with |r| <= ln(2)/2, then exp(r) = (exp(r/2^10))^(2^10)
            const DD ln2{0.6931471805599453, 2.3190468138462996e-17};
            const double k = std::nearbyint(x.hi/ln2.hi);
            DD r = sub(x, mul({k, 0.0}, ln2));
            constexpr int squarings = 10;
            r = ldexp(r, -squarings);
            // expm1(r) by its Taylor series, |r| < 3.4e-4, so 12 terms are more than enough
            DD term = r;
            DD s = r;
            for(int n=2; n<=12; ++n){
                term = div(mul(term, r), {static_cast<double>(n), 0.0});
                s = add(s, term);
            }
            // expm1(2r) = expm1(r)*(expm1(r) + 2): no digits are lost when the result is small
            for(int i=0; i<squarings; ++i){
                s = mul(s, add(s, {2.0, 0.0}));
            }
            return ldexp(add(s, {1.0, 0.0}), static_cast<int>(k));
        }
    }

    // the exact solution of a linear chain 0 -> 1 -> ... -> n-1 with all the n0 atoms in 0 at the start.
    // The decay constants must all be different, the last can be 0 (stable). error[m] is a bound of the rounding
    // error of result[m]: ~1e-29 of the sum of the absolute values of the terms.
    struct Bateman {
        std::vector<double> result;
        std::vector<double> error;
    };

    inline Bateman bateman(const std::vector<double>& lambda, double n0, double t){
        const int n = static_cast<int>(lambda.size());
        Bateman b{std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
        DD production{n0, 0.0};
        for(int m=0; m<n; ++m){
            if(m > 0){
                production = dd::mul(production, {lambda[m - 1], 0.0});
            }
            DD sum{};
            double sum_abs = 0.0;
            for(int j=0; j<=m; ++j){
                DD denominator{1.0, 0.0};
                for(int k=0; k<=m; ++k){
                    if(k != j){
                        // the difference of two doubles is exact in double-double
                        denominator = dd::mul(denominator, dd::two_sum(lambda[k], -lambda[j]));
                    }
                }
                const DD term = dd::div(dd::exp(dd::neg(dd::two_prod(lambda[j], t))), denominator);
                sum = dd::add(sum, term);
                sum_abs += std::abs(term.hi);
            }
            const DD value = dd::mul(production, sum);
            b.result[m] = value.hi + value.lo;
            b.error[m] = 1e-29*std::abs(production.hi)*sum_abs;
        }
        return b;
    }
}

#endif // NIS_TESTS_EXACT_BATEMAN_HPP
