#include <vector>
#include <string>
#include <complex>

namespace NFFT {
using cd = std::complex<double>;
using ld = long double;
using cld = std::complex<ld>;

using ComplexPolynom = std::vector<cld>;

using DoublePolynom = std::vector<ld>;
using IntegerPolynom = std::vector<int>;

const double PI = acos(-1);

cld findWn(size_t n, bool invert) {
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

size_t upperPowerofTwo(size_t n) {
    size_t ans = 1;
    while(ans < n) {
        ans <<= 1;
    }
    return ans;
}

ComplexPolynom toComplex(const IntegerPolynom& p) {
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

ComplexPolynom multiplicate(ComplexPolynom p1, ComplexPolynom p2) {
    size_t n = upperPowerofTwo(p1.size() + p2.size());
    p1.resize(n, cld(0));
    p2.resize(n, cld(0));
    
    fft(p1, false);
    fft(p2, false);
    
    for(size_t i = 0; i < n; ++i) {
        p1[i] *= p2[i];
    }
    
    fft(p1, true);
    
    return p1;
}

IntegerPolynom multiplicate(IntegerPolynom p1, IntegerPolynom p2) {
    ComplexPolynom a1 = toComplex(p1);
    ComplexPolynom a2 = toComplex(p2);
    
    ComplexPolynom a = multiplicate(a1, a2);
    
    return fromComplextoInteger(a);
}

};

std::string block_to_string(int number, int osn) {
    std::string l = std::to_string(number);
    int count_zeros = 0;
    osn /= 10;
    while(number < osn) {
        ++count_zeros;
        osn /= 10;
    }
    if(number == 0) --count_zeros;
    while(count_zeros) {
        l = '0' + l;
        --count_zeros;
    }
    return l;
}

class BigInteger {
    friend bool operator<(const BigInteger&, const BigInteger&);
    friend BigInteger operator-(const BigInteger& number);
    friend std::ostream& operator<<(std::ostream& out, const BigInteger& number);
    friend std::istream& operator>>(std::istream& in, BigInteger& number);
private:
    const int OSN = 10000;
    const int OSN_EXP = 4;
    std::vector<int> blocks;
    bool sign; // true - positive
    
    void normalize_zeros();
    void normalize();
    BigInteger& prefix_increament();
    BigInteger& prefix_decreament();
    BigInteger divide_by_base(size_t order) const;
    BigInteger multiplicate_by_base(size_t order) const;
    
public:
    BigInteger();
    BigInteger(int number);
    BigInteger(const BigInteger& number);
    
    BigInteger& operator=(const BigInteger& number);
    explicit operator bool() const;
    std::string toString() const;
    
    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);
    
    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);
    
};

bool operator<(const BigInteger& first, const BigInteger& second);
bool operator>(const BigInteger& first, const BigInteger& second);
bool operator<=(const BigInteger& first, const BigInteger& second);
bool operator>=(const BigInteger& first, const BigInteger& second);
bool operator==(const BigInteger& first, const BigInteger& second);
bool operator!=(const BigInteger& first, const BigInteger& second);

BigInteger operator-(const BigInteger& number);

BigInteger operator+(const BigInteger& first, const BigInteger& second);
BigInteger operator-(const BigInteger& first, const BigInteger& second);
BigInteger operator*(const BigInteger& first, const BigInteger& second);
BigInteger operator/(const BigInteger& first, const BigInteger& second);
BigInteger operator%(const BigInteger& first, const BigInteger& second);

std::ostream& operator<<(std::ostream& out, const BigInteger& number);
std::istream& operator>>(std::istream& in, BigInteger& number);

BigInteger gcd(const BigInteger& first, const BigInteger& second);
BigInteger abs(const BigInteger& number);
void swap(BigInteger& first, BigInteger& second);


bool operator<(const BigInteger& first, const BigInteger& second) {
    if(first.sign && !second.sign) return false;
    if(!first.sign && second.sign) return true;
    bool delta_differ = false;
    
    // если оба отрицательны, то будем менять все знаки
    if(!first.sign && !second.sign) {
        delta_differ = true;
    }
    
    if(first.blocks.size() != second.blocks.size()) {
        return (first.blocks.size() < second.blocks.size()) != delta_differ;
    }
    size_t i = first.blocks.size();
    int delta;
    while(i > 0) {
        delta = first.blocks[i - 1] - second.blocks[i - 1];
        if(delta || !(i - 1)) return (delta < 0) != delta_differ;
        --i;
    }
    return false;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
       return second < first;
}
   
bool operator==(const BigInteger& first, const BigInteger& second) {
    return !(first < second) && !(second < first);
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return (first < second) || (second < first);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
       return !(first < second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
       return !(second < first);
}

BigInteger operator-(const BigInteger& number) {
    BigInteger result = number;
    if(number != 0) {
        result.sign = !number.sign;
    }
    return result;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result += second;
    return result;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    return first + (-second);
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result *= second;
    return result;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result /= second;
    return result;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result %= second;
    return result;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    char symbol;
    int block = 0;
    int order = 1;
    std::vector<int> queue;
    number.blocks.clear();
    number.sign = true;
    
    in >> std::noskipws;
    while(in >> symbol && !isspace(symbol)) {
        if(symbol == '-') {
            number.sign = false;
            continue;
        }
        queue.push_back(symbol - '0');
    }
    
    for(size_t i = queue.size(); i > 0; --i) {
        block += order * queue[i - 1];
        order *= 10;
        if(order >= number.OSN) {
            number.blocks.push_back(block);
            order = 1;
            block = 0;
        }
    }
    if(order != 1) {
        number.blocks.push_back(block);
    }
    number.normalize_zeros();
    return in;
}


void BigInteger::normalize_zeros() {
    while(blocks.size() > 1 && blocks.back() == 0) {
        blocks.pop_back();
    }
    if(blocks.size() == 1 && blocks.back() == 0) {
        sign = true;
    }
}

void BigInteger::normalize() {
    int next_digit;
    for(size_t i = 0; i < blocks.size(); ++i) {
        if(blocks[i] >= 0) {
            next_digit = blocks[i] / OSN;
            blocks[i] = blocks[i] % OSN;
            if(i == blocks.size() - 1) {
                if(next_digit != 0) {
                    blocks.push_back(next_digit);
                }
                continue;
            }
            blocks[i + 1] += next_digit;
        }
        else{
            if(i == blocks.size() - 1) {
                sign = !sign;
                for(size_t j = 0; j < blocks.size(); ++j) {
                    blocks[j] *= -1;
                }
                normalize();
                return;
            }
            
            next_digit = (-blocks[i]) / OSN;
            blocks[i] = blocks[i] % OSN;
            if(blocks[i] < 0) {
                blocks[i] += OSN;
                ++next_digit;
            }
            blocks[i + 1] -= next_digit;
        }
    }
}

BigInteger& BigInteger::prefix_increament() {
    size_t i = 0;
    ++blocks[i];
    while(blocks[i] >= OSN) {
        blocks[i] -= OSN;
        if(i != blocks.size()) {
            ++blocks[i + 1];
            continue;
        }
        blocks.push_back(1);
    }
    normalize_zeros();
    return *this;
}

BigInteger& BigInteger::prefix_decreament() {
    if(static_cast<bool>(*this)) {
        size_t i = 0;
        --blocks[i];
        while(blocks[i] < 0) {
            blocks[i] += OSN;
            --blocks[i + 1];
        }
    } else {
        ++blocks[0];
        sign = false;
    }
    normalize_zeros();
    return *this;
}

BigInteger BigInteger::divide_by_base(size_t order) const{
    BigInteger res = 0;
    if(order < blocks.size()) res.blocks.clear();
    for(size_t i = order; i < blocks.size(); ++i) {
        res.blocks.push_back(blocks[i]);
    }
    return res;
}

BigInteger BigInteger::multiplicate_by_base(size_t order) const{
    BigInteger res = 0;
    res.blocks.clear();
    for(size_t i = 0; i < order; ++i) {
        res.blocks.push_back(0);
    }
    for(size_t i = 0; i < blocks.size(); ++i) {
        res.blocks.push_back(blocks[i]);
    }
    return res;
}

BigInteger::BigInteger(): blocks({0}), sign(true) {}

BigInteger::BigInteger(int number) {
    sign = true;
    if(number < 0) {
        sign = false;
        number = -number;
    }
    blocks.clear();
    if(number == 0) {
        blocks.push_back(0);
    }
    while(number > 0) {
        blocks.push_back(number % OSN);
        number /= OSN;
    }
}


BigInteger::BigInteger(const BigInteger& number): blocks(number.blocks), sign(number.sign) {}

BigInteger& BigInteger::operator=(const BigInteger& number) {
    blocks = number.blocks;
    sign = number.sign;
    return *this;
}

BigInteger::operator bool() const{
    if(blocks.size() > 1) return true;
    if(blocks[0]) return true;
    return false;
}

std::string BigInteger::toString() const{
    std::string s = "";
    if(!sign) s += '-';
    s += std::to_string(blocks[blocks.size() - 1]);
    for(size_t i = blocks.size() - 1; i > 0; --i) {
        s += block_to_string(blocks[i - 1], OSN);
    }
    return s;
}

BigInteger& BigInteger::operator+=(const BigInteger& number) {
    size_t number_size = number.blocks.size();
    while(blocks.size() < number_size) {
        blocks.push_back(0);
    }
    int delta = (number.sign == sign)? 1 : -1;
    for(size_t i = 0; i < number_size; ++i) {
        blocks[i] += delta * number.blocks[i];
    }
    normalize();
    normalize_zeros();
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& number) {
    return *this += (-number);
}

BigInteger& BigInteger::operator*=(const BigInteger& number) {
    blocks = NFFT::multiplicate(blocks, number.blocks);
    
    sign = (sign == number.sign);
    normalize();
    normalize_zeros();
    return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& number) {
    BigInteger numerator = abs(*this);
    BigInteger denominator = abs(number);
    BigInteger divided;
    BigInteger subtractive;
    
    std::vector<int> digits(0);
    blocks.clear();
    
    int digit;
    int bound_digit;
    int pivot;
    if(numerator.blocks.size() < denominator.blocks.size()){
        *this = 0;
        return *this;
    }
    for(size_t i = numerator.blocks.size() - denominator.blocks.size() + 1; i > 0; --i) {
        divided = numerator.divide_by_base(i - 1);
        digit = 0;
        bound_digit = OSN;
        while(bound_digit > digit + 1) {
            pivot = (digit + bound_digit) / 2;
            ((divided >= pivot * denominator) ? digit : bound_digit) = pivot;
        }
        
        subtractive = digit * denominator;
        numerator -= subtractive.multiplicate_by_base(i - 1);
        digits.push_back(digit);
    }
    
    size_t digit_size = digits.size();
    for(size_t i = 1; i <= digit_size; ++i) {
        blocks.push_back(digits[digit_size - i]);
    }
    if(digit_size == 0){
        blocks.push_back(0);
    }
    
    normalize_zeros();
    sign = (sign == number.sign);

    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& number) {
    BigInteger even = *this / number;
    even *= number;
    *this -= even;
    return *this;
}

BigInteger& BigInteger::operator++() {
    if(!sign) return prefix_decreament();
    return prefix_increament();
}

BigInteger BigInteger::operator++(int) {
    BigInteger temporary = *this;
    ++*this;
    return temporary;
}

BigInteger& BigInteger::operator--() {
    if(!sign) return prefix_increament();
    return prefix_decreament();
}

BigInteger BigInteger::operator--(int) {
    BigInteger temporary = *this;
    --*this;
    return temporary;
}

BigInteger gcd(const BigInteger& first, const BigInteger& second) {
    BigInteger a = abs(first);
    BigInteger b = abs(second);
    if(a == 0 || b == 0) return 0;
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

BigInteger abs(const BigInteger& number) {
    if(number < 0) return -number;
    return number;
}

void swap(BigInteger& first, BigInteger& second) {
    BigInteger temporary = first;
    first = second;
    second = temporary;
}

class Rational {
    friend bool operator<(const Rational& first, const Rational& second);
    friend Rational operator-(const Rational& number);
private:
    BigInteger nominator;
    BigInteger denominator;
    
    void reduction();

public:
    Rational();
    Rational(BigInteger number);
    Rational(int number);
    
    //explicit operator double();
    
    Rational& operator+=(const Rational&);
    Rational& operator-=(const Rational&);
    Rational& operator*=(const Rational&);
    Rational& operator/=(const Rational&);
    
    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;
};

bool operator<(const Rational& first, const Rational& second);
bool operator>(const Rational& first, const Rational& second);
bool operator<=(const Rational& first, const Rational& second);
bool operator>=(const Rational& first, const Rational& second);
bool operator==(const Rational& first, const Rational& second);
bool operator!=(const Rational& first, const Rational& second);

Rational operator-(const Rational& number);

Rational operator+(const Rational& first, const Rational& second);
Rational operator-(const Rational& first, const Rational& second);
Rational operator*(const Rational& first, const Rational& second);
Rational operator/(const Rational& first, const Rational& second);

std::ostream& operator<<(std::ostream& out, const Rational& number) {
    out << number.asDecimal(6);
    return out;
}


Rational::Rational(): nominator(0), denominator(1) {}

Rational::Rational(BigInteger number): nominator(number), denominator(1) {}

Rational::Rational(int number): nominator(number), denominator(1) {}

Rational& Rational::operator+=(const Rational& number) {
    nominator = nominator * number.denominator + denominator * number.nominator;
    denominator = denominator * number.denominator;
    reduction();
    return *this;
}

Rational& Rational::operator-=(const Rational& number) {
    return *this += (-number);
}

Rational& Rational::operator*=(const Rational& number) {
    nominator *= number.nominator;
    denominator *= number.denominator;
    reduction();
    return *this;
}

Rational& Rational::operator/=(const Rational& number) {
    nominator *= number.denominator;
    denominator *= number.nominator;
    if(denominator < 0) {
        denominator = -denominator;
        nominator = -nominator;
    }
    reduction();
    return *this;
}

std::string Rational::toString() const{
    std::string result = "";
    result += nominator.toString();
    if(denominator != 1) {
        result += '/';
        result += denominator.toString();
    }
    return result;
}

std::string Rational::asDecimal(size_t precision) const{
    std::string result = "";
    Rational number = *this;
    BigInteger even;
    bool is_dot = false;
    
    if(number.nominator < 0) {
        result += '-';
        number.nominator = -number.nominator;
    }
    ++precision;
    
    while(precision > 0) {
        even = number.nominator / number.denominator;
        number -= even;
        number *= 10;
        result += even.toString();
        --precision;
        
        if(!is_dot && precision > 0) {
            is_dot = true;
            result += '.';
        }
    }
    return result;
}

bool operator<(const Rational& first, const Rational& second) {
    return first.nominator * second.denominator < second.nominator * first.denominator;
}

bool operator>(const Rational& first, const Rational& second) {
       return second < first;
}
   
bool operator==(const Rational& first, const Rational& second) {
    return !(first < second) && !(second < first);
}

bool operator!=(const Rational& first, const Rational& second) {
    return (first < second) || (second < first);
}

bool operator>=(const Rational& first, const Rational& second) {
       return !(first < second);
}

bool operator<=(const Rational& first, const Rational& second) {
       return !(second < first);
}

Rational operator-(const Rational& number) {
    Rational result = number;
    result.nominator = -result.nominator;
    return result;
}

Rational operator+(const Rational& first, const Rational& second) {
    Rational result = first;
    result += second;
    return result;
}

Rational operator-(const Rational& first, const Rational& second) {
    return first + (-second);
}

Rational operator*(const Rational& first, const Rational& second) {
    Rational result = first;
    result *= second;
    return result;
}

Rational operator/(const Rational& first, const Rational& second) {
    Rational result = first;
    result /= second;
    return result;
}

void Rational::reduction() {
    BigInteger coefficient = gcd(nominator, denominator);
    if(coefficient == 0){
        denominator = 1;
        return;
    }
    nominator /= coefficient;
    denominator /= coefficient;
}

