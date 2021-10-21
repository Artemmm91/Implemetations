#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>

class Complex {
private:
    long double r, im;
public:
    Complex(): r(0), im(0) {}
    Complex(long double real, long double im): r(real), im(im) {}
    Complex(long double real): r(real), im(0) {}
    Complex(const Complex& another) = default;

    long double real() const{
        return r;
    }

    long double imag() const {
        return im;
    }

    Complex& operator+=(const Complex& another){
        r += another.r;
        im += another.im;
        return *this;
    }

    Complex& operator*=(const Complex& another){
        Complex it(r * another.r - im * another.im, r * another.im + im * another.r);
        *this = it;
        return *this;
    }

    Complex operator+(const Complex& another) const{
        return Complex(r + another.r, im + another.im);
    }

    Complex operator*(const Complex& another) const{
        return Complex(r * another.r - im * another.im, r * another.im + im * another.r);
    }

    Complex operator-() const{
        return Complex(-r, -im);
    }

    Complex operator-(const Complex& another) const{
        return Complex(r - another.r, im - another.im);
    }

    Complex& operator/=(size_t n){
        r /= n;
        im /= n;
        return *this;
    }

    Complex conj() const{
        return Complex(r, -im);
    }
};

void swap(Complex& this_c, Complex& an_c){
    Complex tmp = this_c;
    this_c = an_c;
    an_c = tmp;
}

namespace NFFT {
    using ld = long double;
    using cld = Complex;
    using ComplexPolynom = std::vector<cld>;
    using DoublePolynom = std::vector<long double>;
    using IntegerPolynom = std::vector<long long>;
    double PI = 2 * acos(-1);

    size_t reverse(size_t num, size_t lg_n) {
        size_t res = 0;
        for (size_t i = 0; i < lg_n; i++) {
            if (num & (1 << i))
                res |= 1 << (lg_n - 1 - i);
        }
        return res;
    }

    void revSort(ComplexPolynom& a){
        size_t n = a.size();
        size_t log2n = log2(n);

        for (size_t i = 0; i < n; ++i) {
            size_t rev_i = reverse(i, log2n);
            if (i < rev_i)
                swap(a[i], a[rev_i]);
        }
    }

    void fftNonRecursive(ComplexPolynom& a, bool invert) {
        size_t n = a.size();
        int delta = invert ? -1 : 1;
        PI *= delta;

        revSort(a);

        for (size_t len = 2; len <= n; len <<= 1) {
            double ang = PI / len;
            cld wlen(cos(ang), sin(ang));
            for (int i = 0; i < n; i += len) {
                cld w(1);
                for (int j = 0; j < len / 2; ++j) {
                    cld u = a[i+j];
                    cld v = a[i+j+len/2] * w;
                    a[i+j] = u + v;
                    a[i+j+len/2] = u - v;
                    w *= wlen;
                }
            }
        }

        if (invert) {
            for (cld& x : a)
                x /= n;
        }
    }

    size_t upperPowerofTwo(size_t n) {
        size_t ans = 1;
        while (ans < n) {
            ans <<= 1;
        }
        return ans;
    }

    ComplexPolynom multiplicate(const ComplexPolynom& p1, const ComplexPolynom& p2, size_t n) {
        ComplexPolynom D(n), Ans(n);
        for(size_t i = 0; i < n; ++i){
            D[i] = cld(p1[i].real(), p2[i].real());
        }

        fftNonRecursive(D, false);

        Ans[0] = D[0].real() * D[0].imag();

        for (size_t i = 1; i < n; ++i) {
            cld a_w = (D[i] + D[n - i].conj()) * cld(0.5);
            cld b_w = (D[i] - D[n - i].conj()) * cld(0, -0.5);
            Ans[i] = a_w * b_w;
        }

        fftNonRecursive(Ans, true);

        return Ans;
    }

    IntegerPolynom multiplicate(const IntegerPolynom& a1, const IntegerPolynom& a2){
        ComplexPolynom p1(a1.begin(), a1.end());
        ComplexPolynom p2(a2.begin(), a2.end());
        size_t n = upperPowerofTwo(p1.size() + p2.size());
        p1.resize(n, cld(0));
        p2.resize(n, cld(0));
        
        ComplexPolynom ans = multiplicate(p1, p2, n);
        
        IntegerPolynom p(n);
        for (size_t i = 0; i < n; ++i) {
            p[i] = std::round(ans[i].real());
        }
        return p;
    }

    DoublePolynom multiplicate(const DoublePolynom& a1, const DoublePolynom& a2){
        ComplexPolynom p1(a1.begin(), a1.end());
        ComplexPolynom p2(a2.begin(), a2.end());
        size_t n = upperPowerofTwo(p1.size() + p2.size());
        p1.resize(n, cld(0));
        p2.resize(n, cld(0));
        
        ComplexPolynom ans = multiplicate(p1, p2, n);
        
        DoublePolynom p(n);
        for (size_t i = 0; i < n; ++i) {
            p[i] = ans[i].real();
        }
        return p;
    }

    DoublePolynom difference(const DoublePolynom* a, const DoublePolynom* b){
        DoublePolynom c(std::max(a->size(), b->size()));
        if(a->size() > b->size()){
            swap(a, b);
        }
        for(size_t i = 0; i < a->size(); ++i){
            c[i] = (*a)[i] - (*b)[i];
        }
        for(size_t i = a->size(); i < b->size(); ++i){
            c[i] = -(*b)[i];
        }
        return c;
    }

    DoublePolynom reciprocal(const DoublePolynom& a, long long precision){
        if(precision == 1) return {1 / a[0]};
        DoublePolynom bk = reciprocal(a, precision / 2);
        DoublePolynom a2 = {2};
        DoublePolynom abk = multiplicate(bk, a);
        DoublePolynom answer = difference(&a2, &abk);
        answer = multiplicate(answer, bk);
        answer.resize(precision);
        return answer;
    }
    
    size_t normal_size(const DoublePolynom& a){
        size_t n = a.size();
        size_t i = n;
        while(i > 0 && a[i - 1] == 0){
            --i;
        }
        return i;
    }

    DoublePolynom divide(const DoublePolynom& a, const DoublePolynom& b){
        DoublePolynom Ar(a.begin(), a.end()), Br(b.begin(), b.end());
        size_t n = normal_size(Ar), m = normal_size(Br);
        Ar.resize(n);
        Br.resize(m);
        std::reverse(Ar.begin(), Ar.end());
        std::reverse(Br.begin(), Br.end());
        DoublePolynom br_1 = reciprocal(Br, upperPowerofTwo(n - m + 1));
        br_1.resize(n - m + 1);
        DoublePolynom q = multiplicate(Ar, br_1);
        q.resize(n - m + 1);
        std::reverse(q.begin(), q.end());
        return q;
    }

    DoublePolynom remainder(const DoublePolynom& a, const DoublePolynom& b){
        IntegerPolynom q = divide(a, b);
        IntegerPolynom qb = multiplicate(q, b);
        IntegerPolynom r = difference(&a, &qb);
        r.resize(normal_size(r));
        return r;
    }

}; // namespace NFFT
