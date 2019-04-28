//
//  uint_t.h
//  uint_t
//
//  Created by Arthur Sun on 4/26/19.
//  Copyright © 2019 Arthur Sun. All rights reserved.
//

#ifndef uint_t_h
#define uint_t_h

#ifndef uint_t_initial_capacity
#define uint_t_initial_capacity 8
#endif

#include <stack>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

class uint_t
{
    
public:
    
    // chunk type
    using uintk = unsigned int;
    
    // (2 * width) of 'uintk'
    using uintk2 = unsigned long;
    
    static const uintk mask = ~uintk(0);
    
    // must be a power of 2
    static const uintk bits = 32;
    
    static const uintk2 base = (uintk2(1) << bits);
    
    static const uintk karatsuba_threshold = 32;
    
    static std::vector<uint_t> powerCache[1 << CHAR_BIT];
    
private:
    
    uintk *data;
    
    uintk size;
    
    uintk capacity;
    
    /// returns radix ^ (2 ^ exp)
    static const uint_t& getPower(uintk radix, uintk exp) {
        std::vector<uint_t>& line = powerCache[radix];
        
        if(line.size() == 0)
            line.push_back(radix);
        
        for(uintk i = (uintk)line.size(); i <= exp; ++i) {
            line[i] = (line[i - 1] << 1);
        }
        
        return line[exp];
    }
    
    /// exclusive
    inline uint_t lowerBits(uintk position) const {
        return uint_t(data, data + position);
    }
    
    /// inclusive
    inline uint_t upperBits(uintk position) const {
        return uint_t(data + position, data + size);
    }
    
    // left shifts the whole number (bits * bbits) bits
    void shiftLeftBits(uintk bbits) {
        if(isUint(0)) return;
        if(bbits == 0) return;
        
        realloc(bbits + size);
        
        for(uintk i = 0; i < size; ++i) {
            data[size - i - 1 + bbits] = data[size - i - 1];
            data[size - i - 1] = 0;
        }
        
        size += bbits;
    }
    
    void shiftRightBits(uintk bbits) {
        if(isUint(0)) return;
        if(bbits == 0) return;
        
        if(bbits >= size) {
            size = 0;
            return;
        }
        
        size -= bbits;
        
        for(uintk i = 0; i < size; ++i) {
            data[i] = data[i + bbits];
        }
    }
    
    // for bbits < bits
    void shiftLeftSingle(uintk bbits) {
        if(isUint(0)) return;
        if(bbits == 0) return;
        
        uintk left = 0;
        
        for(uintk i = 0; i < size; ++i) {
            uintk _left = data[i] >> (bits - bbits);
            data[i] <<= bbits;
            data[i] |= left;
            left = _left;
        }
        
        if(left > 0) {
            realloc(size + 2);
            data[size] = left;
            ++size;
        }
    }
    
    void shiftRightSingle(uintk bbits) {
        if(isUint(0)) return;
        if(bbits == 0) return;
        
        if(size == 1) {
            data[0] >>= bbits;
            if(data[0] == 0)
                size = 0;
            return;
        }
        
        uintk right = 0;
        
        for(uintk i = 0; i < size; ++i) {
            uintk _right = data[size - i - 1] << (bits - bbits);
            data[size - i - 1] >>= bbits;
            data[size - i - 1] |= right;
            
            right = _right;
            
            if(i == 0 && data[size - i - 1] == 0)
                --size;
        }
    }
    
    std::string toString(uintk i, uintk radix) const {
        std::string result;
        uintk n = data[i];
        while(n > 0) {
            result.insert(result.begin(), (n % radix) + '0');
            n /= radix;
        }
        return result;
    }
    
    static void toString(const uint_t& x, uintk radix, int digits, std::string* str) {
        if(x.getSize() == 0)
            return;
        
        if(x.getSize() < 2) {
            std::string s = x.toString(0, radix);
            
            if ((s.length() < digits) && (str->length() > 0)) {
                for (size_t i = s.length(); i < digits; i++) {
                    str->push_back('0');
                }
            }
            
            str->append(s);
            
            return;
        }
        
        uintk bitLength = x.size * bits;
        
        uintk _n = log(bitLength * M_LN2 / log(radix)) / M_LN2 - 0.5f;
        
        uint_t p = getPower(radix, _n);
        /*
        auto pair = x.quotientAndReminder(p);
        
        int k = 1 << _n;
        
        toString(pair.first, radix, digits - k, str);
        toString(pair.second, radix, k, str);
         */
    }

public:
    
    uint_t() : size(0), capacity(uint_t_initial_capacity) {
        data = (uintk*)malloc(sizeof(uintk) * capacity);
    }
    
    uint_t(const uint_t& x) : size(x.getSize()), capacity(x.getCapacity()) {
        data = (uintk*)malloc(sizeof(uintk) * capacity);
        memcpy(data, x.begin(), sizeof(uintk) * size);
    }
    
    uint_t(const uintk* begin, const uintk* end) {
        size = (uintk)(end - begin);
        capacity = size * 2 + uint_t_initial_capacity;
        data = (uintk*)malloc(sizeof(uintk) * capacity);
        memcpy(data, begin, sizeof(uintk) * size);
    }
    
    template <typename T>
    uint_t(const T& x) : size(0), capacity(uint_t_initial_capacity) {
        data = (uintk*)malloc(sizeof(uintk) * capacity);
        uintk i = 0;
        while (true) {
            if((sizeof(T) * CHAR_BIT) <= i || (x >> i) == 0)
                break;
            push((x >> i) & mask);
            i += bits;
        }
    }
    
    ~uint_t() {
        free(data);
    }
    
    uint_t& operator = (const uint_t& x) {
        size = x.getSize();
        realloc(x.getCapacity());
        memcpy(data, x.begin(), sizeof(uintk) * size);
        return *this;
    }
    
    void realloc(uintk _cap) {
        if(_cap > capacity) {
            capacity = _cap;
            data = (uintk*)::realloc(data, sizeof(uintk) * capacity);
        }
    }
    
    void push(uintk x) {
        if(size >= capacity) {
            realloc(capacity * 2);
        }
        data[size++] = x;
    }
    
    uint_t& shiftLeft(uintk _bits) {
        shiftLeftBits(_bits / bits);
        shiftLeftSingle(_bits & (bits - 1));
        return *this;
    }
    
    uint_t& shiftRight(uintk _bits) {
        shiftRightBits(_bits / bits);
        shiftRightSingle(_bits & (bits - 1));
        return *this;
    }
    
    std::string toString(uintk radix) const {
        if(radix >= (1 << CHAR_BIT))
            radix = 10;
        
        std::string str;
        toString(*this, radix, 0, &str);
        return str;
    }
    
    inline uintk getSize() const {
        return size;
    }
    
    inline uintk getCapacity() const {
        return capacity;
    }
    
    inline const uintk* begin() const {
        return data;
    }
    
    inline const uintk* end() const {
        return data + size;
    }
    
    inline uintk*& begin() {
        return data;
    }
    
    inline uintk operator [] (uintk i) const {
        return data[i];
    }
    
    inline uintk& operator [] (uintk i) {
        return data[i];
    }
    
    inline bool isUint (uintk x) const {
        return size == 0 || (size == 1 && data[0] == x);
    }
    
    inline uintk uint () const {
        if(size == 0)
            return 0;
        
        return data[0];
    }
    
    inline uintk2 uint2 () const {
        if(size == 0)
            return 0;
        
        if(size == 1)
            return data[0];
        
        return data[0] + data[1] * base;
    }
    
    inline bool isEven () const {
        return size == 0 || (data[0] & 1) == 0;
    }
    
    inline bool isOdd () const {
        return !isEven();
    }
    
    
    void multiplyByUint (uintk b) {
        if(size == 0 || b == 0) {
            size = 0;
            return;
        }
        
        uintk carry = 0;
        for(uintk i = 0; i < size; ++i) {
            uintk2 product = uintk2(data[i]) * b + carry;
            data[i] = product & mask;
            carry = product >> bits;
        }
        
        if(carry > 0)
            push(carry);
    }
    
    void add (const uint_t& b) {
        uintk carry = 0;
        uintk i = 0;
        
        uintk sizeB = b.getSize();
        
        uintk2 sum;
        
        while(true) {
            if(i >= size && i >= sizeB) {
                if(carry != 0)
                    data[i] = carry;
                
                return;
            }
            
            if(i >= size) {
                if(carry == 0) {
                    realloc(sizeB);
                    size = sizeB;
                    memcpy(data + i, b.begin() + i, sizeof(uintk) * (sizeB - i));
                    return;
                }
                sum = uintk2(b[i]) + carry;
            }else if(i >= sizeB) {
                if(carry == 0)
                    return;
                
                sum = uintk2(data[i]) + carry;
            }else{
                sum = uintk2(data[i]) + b[i] + carry;
            }
            
            data[i] = (sum & mask);
            carry = sum >> bits;
            ++i;
        }
        
        size = i;
    }
    
    void subtract (const uint_t& b) {
        uintk borrow = 0;
        uintk i = 0;
        
        uintk sizeB = b.getSize();
        
        if(size == 0 || size < sizeB) {
            size = 0;
            return;
        }
        
        if(b.isUint(0))
            return;
        
        while(true) {
            if(i >= sizeB && borrow == 0) {
                if(i >= size)
                    size = i - 1;
                return;
            }
            
            uintk2 td = uintk2(i >= sizeB ? 0 : b[i]) + borrow;
            
            if(data[i] < td) {
                if(i >= size - 1) {
                    size = 0;
                    return;
                }
                
                data[i] = (uintk)(base + data[i] - td);
                borrow = 1;
            }else{
                data[i] = (uintk)(data[i] - td);
                borrow = 0;
            }
            
            ++i;
        }
        
        size = i;
    }
    
    
    
    
    /// operators
    
    inline uint_t& operator += (const uint_t& a) {
        add(a);
        return *this;
    }
    
    inline uint_t& operator -= (const uint_t& a) {
        subtract(a);
        return *this;
    }
    
    inline uint_t& operator ++ () {
        add(1);
        return *this;
    }
    
    inline uint_t& operator -- () {
        subtract(1);
        return *this;
    }
    
    inline uint_t& operator *= (uintk a) {
        multiplyByUint(a);
        return *this;
    }
    
    inline uint_t& operator *= (const uint_t& a) {
        return (*this = *this * a);
    }
    
    inline uint_t operator << (uintk x) {
        return shiftLeft(x);
    }
    
    inline uint_t operator >> (uintk x) {
        return shiftRight(x);
    }
    
    /**
     * dividing digit by digit
     */
    std::pair<uint_t, uintk> quotientAndReminderUint (uintk b) const {
        std::pair<uint_t, uintk> result;
        
        uintk mod = 0;
        uintk i = 0;
        
        int startI = -1;
        
        while(i < size) {
            uintk2 mb = uintk2(mod) << bits;
            
            uintk2 k = (uintk2(data[i]) + mb) / b;
            
            if(k > 0 && startI == -1)
                startI = (int)i;
            
            if(startI != -1)
                result.first.push((uintk)k);
            
            mod = (mb + data[i]) % b;
            
            ++i;
        }
        
        result.second = mod;

        return result;
    }
    
    std::pair<uint_t, uint_t> quotientAndReminderN (const uint_t& b) const {
        std::pair<uint_t, uint_t> result;
        
        uintk m = 1l << (bits - 1);
        
        return result;
    }
    
    std::pair<uint_t, uint_t> quotientAndReminder (const uint_t& b) const {
        char comp = compare(*this, b);
        
        if(b.isUint(0)) {
            throw "divisor is zero";
            return std::pair<uint_t, uint_t>(0, 0);
        }
        
        if(comp == 0)
            return std::pair<uint_t, uint_t>(1, 0);
        
        if(comp == -1)
            return std::pair<uint_t, uint_t>(0, 0);
        
        std::pair<uint_t, uint_t> result;
        
        
        
        return result;
    }
    
    
    /// friend functions
    
    
    /// return sign(a - b)
    friend char compare (const uint_t& a, const uint_t& b) {
        if(a.getSize() > b.getSize())
            return 1;
        
        if(a.getSize() < b.getSize())
            return -1;
        
        uintk n = a.getSize();
        
        for(uintk i = 0; i < n; ++i) {
            if(a[n - i - 1] > a[n - i - 1])
                return 1;
            
            if(a[n - i - 1] < a[n - i - 1])
                return -1;
        }
        
        return 0;
    }
    
    
    
    /**
     * adding each digit in base (2 ^ bits)
     */
    friend uint_t operator + (const uint_t& a, const uint_t& b) {
        uint_t result(a);
        result.add(b);
        return result;
    }
    
    /**
     * subtracting each digit in base (2 ^ bits)
    */
    friend uint_t operator - (const uint_t& a, const uint_t& b) {
        uint_t result(a);
        result.subtract(b);
        return result;
    }
    
    friend uint_t multiplyUint (const uint_t& a, uintk b) {
        uint_t result = a;
        result.multiplyByUint(b);
        return result;
    }
    
    
    friend uint_t multiply (const uint_t& a, const uint_t& b) {
        uint_t result;
        
        uintk an = a.getSize();
        uintk bn = b.getSize();
        
        if(an == 0 || bn == 0)
            return uint_t();
        
        uintk size = an + bn;
        
        result.realloc(size);
        
        result.size = size - 1;
        
        uintk carry = 0;
        for (uintk i = 0; i < bn; ++i) {
            uintk2 product = uintk2(a[0]) * b[i] + carry;
            result[i] = product & mask;
            carry = product >> bits;
        }
        
        result[size - an] = carry;
        
        for (uintk i = 1; i < an; ++i) {
            carry = 0;
            for (uintk j = 0; j < bn; ++j) {
                uintk k = i + j + 1;
                uintk2 product = uintk2(a[i]) * b[j] + result[k] + carry;
                result[k] = product & mask;
                carry = product >> bits;
            }
            result[i] = carry;
        }
        
        if(result[size - 1] != 0)
            ++result.size;
        
        return result;
    }
    
    
    friend uint_t karatsubaMultiply (const uint_t& a, const uint_t& b) {
        uintk half = (std::max(a.getSize(), b.getSize()) + 1) / 2;
        
        uint_t lowerA = a.lowerBits(half);
        uint_t upperA = a.upperBits(half);
        uint_t lowerB = b.lowerBits(half);
        uint_t upperB = b.upperBits(half);
        
        uint_t p0 = upperA * upperB;
        uint_t p1 = lowerA * lowerB;
        uint_t p2 = (lowerA + upperA) * (lowerB + upperB);
        
        return (p0 << (bits * 2 * half)) + ((p2 - p1 - p0) << (bits * half)) + p1;
    }
    
    friend uint_t operator * (const uint_t& a, const uint_t& b) {
        if(a.getSize() == 0 || b.getSize() == 0)
            return 0;
        
        if(a.getSize() == 1)
            return multiplyUint(b, a.uint());
        
        if(b.getSize() == 1)
            return multiplyUint(a, b.uint());
        
        if(std::min(a.getSize(), b.getSize()) < karatsuba_threshold)
            return multiply(a, b);
        
        return karatsubaMultiply(a, b);
    }
};

std::vector<uint_t> uint_t::powerCache[1 << CHAR_BIT];

#endif /* uint_t_h */
