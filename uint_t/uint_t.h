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
    
    static const uintk2 base = uintk2(1) << bits;
    
private:

    uintk *data;
    
    uintk size;
    
    uintk capacity;
    
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
            left = data[i] >> (bits - bbits);
            data[i] <<= bbits;
            data[i] |= left;
        }
        
        if(left > 0) {
            realloc(size * 2 + 1);
            data[size] = left;
            ++size;
        }
    }
    
    void shiftRightSingle(uintk bbits) {
        if(isUint(0)) return;
        if(bbits == 0) return;
        
        if(size == 1) {
            data[0] >>= bbits;
            return;
        }
        
        uintk right = 0;
        
        for(uintk i = 0; i < size; ++i) {
            right = data[size - i - 1] << (bits - bbits);
            data[size - i - 1] >>= bbits;
            data[size - i - 1] |= right;
            
            if(i == 0 && data[size - i - 1] == 0)
                --size;
        }
    }
    
public:
    
    uint_t() : size(0), capacity(uint_t_initial_capacity) {
        data = (uintk*)malloc(sizeof(uintk) * capacity);
    }
    
    uint_t(const uint_t& x) : size(x.getSize()), capacity(x.getCapacity()) {
        data = (uintk*)malloc(sizeof(uintk) * capacity);
        memcpy(data, x.begin(), sizeof(uintk) * size);
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
    
    std::string toString() const {
        return "0";
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
    
    inline bool isUint (uintk x) const {
        return size == 0 || (size == 1 && data[0] == x);
    }
    
    inline uintk uint () const {
        if(size == 0)
            return 0;
        
        return data[0];
    }
    
    inline bool isEven () const {
        return size == 0 || (data[0] & 1) == 0;
    }
    
    inline bool isOdd () const {
        return !isEven();
    }
    
    
    
    
    /// operators
    
    inline uint_t& operator += (const uint_t& a) {
        return (*this = *this + a);
    }
    
    inline uint_t& operator -= (const uint_t& a) {
        return (*this = *this - a);
    }
    
    inline uint_t& operator ++ () {
        return (*this = *this + 1);
    }
    
    inline uint_t& operator -- () {
        return (*this = *this - 1);
    }
    
    inline uint_t& operator *= (const uint_t& a) {
        return (*this = *this * a);
    }
    
    inline uint_t operator << (uintk x) const {
        return uint_t(*this).shiftLeft(x);
    }
    
    inline uint_t operator >> (uintk x) const {
        return uint_t(*this).shiftRight(x);
    }
    
    
    
    
    /// friend functions
    
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
     * adding each digit in base (2 ^ 32)
     */
    friend uint_t operator + (const uint_t& a, const uint_t& b) {
        uint_t result;
        
        uintk carry = 0;
        uintk i = 0;
        
        uintk sizeA = a.getSize();
        uintk sizeB = b.getSize();
        
        uintk2 sum;
        
        while(true) {
            if(i >= sizeA && i >= sizeB) {
                if(carry != 0) result.push(carry);
                break;
            }
            
            if(i >= sizeA) {
                if(carry == 0) {
                    result.realloc(b.getSize());
                    result.size = b.getSize();
                    memcpy(result.begin() + i, b.begin() + i, sizeof(uintk) * (b.getSize() - i));
                    break;
                }
                sum = uintk2(b[i]) + carry;
            }else if(i >= sizeB) {
                if(carry == 0) {
                    result.realloc(a.getSize());
                    result.size = a.getSize();
                    memcpy(result.begin() + i, a.begin() + i, sizeof(uintk) * (a.getSize() - i));
                    break;
                }
                sum = uintk2(a[i]) + carry;
            }else{
                sum = uintk2(a[i]) + b[i] + carry;
            }

            result.push(sum & mask);
            carry = sum > mask;
            ++i;
        }
        
        return result;
    }
    
    /**
     * subtracting each digit in base (2 ^ 32)
    */
    friend uint_t operator - (const uint_t& a, const uint_t& b) {
        uint_t result;
        
        uintk borrow = 0;
        uintk i = 0;
        
        uintk sizeA = a.getSize();
        uintk sizeB = b.getSize();
        
        if(sizeA < sizeB)
            return result;
        
        while(true) {
            if(i >= sizeB && borrow == 0) {
                if(i >= sizeA) break;
                result.realloc(a.getSize());
                result.size = a.getSize();
                memcpy(result.begin() + i, a.begin() + i, sizeof(uintk) * (a.getSize() - i));
                break;
            }
            
            uintk2 td = uintk2(i >= sizeB ? 0 : b[i]) + borrow;
            
            if(a[i] < td) {
                if(i >= sizeA - 1)
                    return uint_t();
                
                result.push((uintk)(base + a[i] - td));
                borrow = 1;
            }else{
                result.push((uintk)(a[i] - td));
                borrow = 0;
            }
            
            ++i;
        }
        
        return result;
    }
    
    
    /**
     * if number is odd, minus 1
     * else divide by 2
     * repeat
     
     * apply this to number 'b'
     * 123 -> 122 -> 61 -> 60 -> 30 -> 15 -> 14 -> 7 -> 6 -> 3 -> 2 -> 1
     *     +1 <-  *2 <- +1   <-  *4   <-  +1 <- *2 <- ...
     
     * save the steps
     
     * apply the steps to number 'a'
     
     * now you have 'a' times 'b'
     */
    friend uint_t operator * (const uint_t& a, const uint_t& b) {
        // base cases
        if(a.isUint(0) || b.isUint(0)) return uint_t();
        
        if(b.isUint(1)) return a;
        
        if(b.isUint(2)) return (a << 1);
        
        if(a.isUint(1)) return b;
        
        if(a.isUint(2)) return (b << 1);
        
        std::stack<uintk> work;
        
        uint_t _b(b);
        
        while(!_b.isUint(1)) {
            if(_b.isEven()) {
                if(work.empty() || (work.top() == 0)) {
                    work.push(1);
                }else{
                    ++work.top();
                }
                _b.shiftRight(1);
            }else{
                --_b;
                work.push(0);
            }
        }
        
        uint_t result(a);
        
        while(!work.empty()) {
            uintk q = work.top();
            if(q == 0) {
                result += a;
            }else{
                result.shiftLeft(q);
            }
            work.pop();
        }
        
        return result;
    }
};

#endif /* uint_t_h */
