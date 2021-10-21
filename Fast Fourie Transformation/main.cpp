#include <iostream>
#include <complex>
#include <vector>

namespace NFFT{
using cd = std::complex<double>;
using ld = long double;
using cld = std::complex<ld>;

using ComplexPolynom = std::vector<cld>;

using DoublePolynom = std::vector<ld>;
using IntegerPolynom = std::vector<int>;

const double PI = acos(-1);

cld findWn(size_t n, bool invert){
    ld ang = 2 * PI / n * (invert ? -1 : 1);
    return cld(cosl(ang), sinl(ang));
}

void fft(ComplexPolynom& a, bool invert) {
    size_t n = a.size();
    
    if (n == 1) return;

    ComplexPolynom a0(n / 2), a1(n / 2);
    for (size_t i = 0; 2 * i < n; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }
    
    fft(a0, invert);
    fft(a1, invert);

    cld w(1);
    cld wn = findWn(n, invert);
    
    for (size_t i = 0; 2 * i < n; ++i) {
        cld u = a0[i];
        cld v = w * a1[i];
        
        a[i] = u + v;
        a[i + n/2] = u - v;
        
        // в итоге после всех рекурсий поделится на n
        if (invert) {
            a[i] /= 2;
            a[i + n/2] /= 2;
        }
        
        w *= wn;
    }
}

size_t upperPowerofTwo(size_t n){
    size_t ans = 1;
    while(ans < n){
        ans <<= 1;
    }
    return ans;
}

ComplexPolynom toComplex(const IntegerPolynom& p){
    ComplexPolynom a(p.size());
    for(size_t i = 0; i < p.size(); ++i) {
        a[i] = cld(p[i]);
    }
    return a;
}

IntegerPolynom fromComplextoInteger(const ComplexPolynom& a) {
    IntegerPolynom p(a.size());
    for(size_t i = 0; i < a.size(); ++i) {
        p[i] = std::round(a[i].real());
    }
    return p;
}

ComplexPolynom multiplicate(ComplexPolynom p1, ComplexPolynom p2){
    size_t n = upperPowerofTwo(p1.size() + p2.size());
    p1.resize(n, cld(0));
    p2.resize(n, cld(0));
    
    fft(p1, false);
    fft(p2, false);
    
    for(size_t i = 0; i < n; ++i){
        p1[i] *= p2[i];
    }
    
    fft(p1, true);
    
    return p1;
}

IntegerPolynom multiplicate(IntegerPolynom p1, IntegerPolynom p2){
    ComplexPolynom a1 = toComplex(p1);
    ComplexPolynom a2 = toComplex(p2);
    
    ComplexPolynom a = multiplicate(a1, a2);
    
    return fromComplextoInteger(a);
}

};

