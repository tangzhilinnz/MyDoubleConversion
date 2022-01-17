#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

extern void nasonov9(char* p, uint32_t u);

typedef union {
  double n;
  uint64_t u64;
  struct {
    uint32_t lo;
    uint32_t hi;
  } u32;
} TValue;

void decode(double n) {
    TValue t;
    t.n = n;
    if ((t.u32.hi << 1) >= 0xffe00000) {
        // binary form: 
        // t.u32.hi = Bx-111-1111-1111-xxxx-xxxx-xxxx-xxxx-xxxx
        // t.u32.hi << 1 = B1111-1111-111x-xxxx-xxxx-xxxx-xxxx-xxxx
        if (((t.u32.hi & 0x000fffff) | t.u32.lo) != 0) {
            printf("NaN\n");
        }
        else {
            printf("Infinity\n");
        }
    }
    else {
        int32_t e = (t.u32.hi >> 20) & 0x7ff;
        uint64_t m = t.u32.hi & 0xfffff;
        if (e == 0) {
            e++;
        }
        else {
            m |= 0x100000;
        }
        // if (t.u32.lo == 0)
        //     n = (m * 2^32 + t.u32.lo) * 2^(e - 1075) 
        //       = m * 2^32 * 2^(e - 1075) = m * 2^(e - 1043)
        // else 
        //     n = (m * 2^32 + t.u32.lo) * 2^(e - 1075)
        e -= 1043;
        if (t.u32.lo) {
            e -= 32;
            m = (m << 32) | t.u32.lo;
        }
        printf("%llu * 2^%d\n", (long long unsigned)m, (int)e);
    }
}

void nd_print(char* p, uint32_t* nd, int32_t ndlo, int32_t ndhi) {
    int32_t i;
    for (i = ndhi; i >= 0; --i) {
        nasonov9(p, nd[i & 127]); p += 9;
    }
    *p++ = '.';
    for (; i >= ndlo; --i) {
        nasonov9(p, nd[i & 127]); p += 9;
    }
    *p = 0;
}

// Multiplying by 2e takes a bit of effort. Let's start by considering the case
// where e is negative, at which point we're really dividing by 2^-e. In turn, 
// this is dividing by 2 for -e times. Just dividing by 2 first:
int32_t nd_div2(uint32_t* nd, int32_t ndlo, int32_t ndhi) {
    uint32_t i = ndhi & 127;
    uint32_t carry = 0;

    for (;;) {
        uint32_t val = nd[i];

        // (val & 1) is the first bit of the binary form of piece_i before 
        // dividing by 2, 
        // (nd[i] * (10^9)^i) / 2 
        //   = (val - (val & 1) + (val & 1)) / 2 * (10^9)^i
        //   = (val >> 1) * (10^9)^i + (val & 1) / 2 * (10^9)^i
        //   = (val >> 1) * (10^9)^i + (val & 1) * 500000000 * (10^9)^(i-1)
        nd[i] = (val >> 1) + carry;
        carry = (val & 1) * 500000000;
        if (i == (ndlo & 127)) break;
        i = (i - 1) & 127;
    }

    if (carry) nd[--ndlo & 127] = carry;
    return ndlo;
}

// We can generalise this to dividing by 2^k:
int32_t nd_div2k_9(uint32_t* nd, int32_t ndlo, int32_t ndhi, uint32_t k) {
    uint32_t mask = (1U << k) - 1;
    uint32_t mul = 1000000000 >> k;
    uint32_t i = ndhi & 127;
    uint32_t carry = 0;

    for (;;) {
        uint32_t val = nd[i];

        // (val & mask) is the first k bits of the binary form of piece_i before 
        // dividing by 2^k, 
        // (nd[i] * (10^9)^i) / 2^k 
        //   = (val - (val & mask) + (val & mask)) / 2^k * (10^9)^i
        //   = (val >> k) * (10^9)^i + (val & mask) / 2^k * (10^9)^i
        //   = (val >> k) * (10^9)^i + (val & mask) * ((10^9) >> k) * (10^9)^(i-1)
        nd[i] = (val >> k) + carry;
        carry = (val & mask) * mul;
        // carry_max 
        //   = mask * mul 
        //   = (2^k - 1) * (10^9 / 2^k) 
        //   = 10^9 - 10^9 / 2^k < 10^9
        // (val_max >> k) + carry_max
        //   = (10^9 - 1) / 2^k + (2^k - 1) * (10^9 / 2^k)
        // k = 1, (val_max >> 1) + carry_max 
        //          = 499_999_999 + 500_000_000 = 999_999_999
        // k = 2, (val_max >> 2) + carry_max 
        //          = 249_999_999 + 750_000_000 = 999_999_999
        // k = 3, (val_max >> 3) + carry_max 
        //          = 124_999_999 + 875_000_000 = 999_999_999
        // ... ... ... ... ... ... ... ... ... ... ... ... ...
        // k = 9, (val_max >> 9) + carry_max 
        //          = 1_953_124 + 998_046_875 = 999_999_999
        if (i == (ndlo & 127)) break;
        i = (i - 1) & 127;
    }

    if (carry) nd[--ndlo & 127] = carry;
    return ndlo;
}

// The above will work for k between zero and nine inclusive; if k is larger 
// than 9, then 1000000000 >> k can no longer be represented by an integer.
// (1000000000 / 2^9 = 1953125) 1953125 cannot be divided by 2, so 1000000000
// cannot be divided by 2^k (k > 9).
// As such, the complete solution has to start by dividing in batches of 2^9:
int32_t nd_div2k(uint32_t* nd, int32_t ndlo, int32_t ndhi, uint32_t k) {
    while (k >= 9) {
        uint32_t i = ndhi & 127;
        uint32_t carry = 0;

        for (;;) {
            uint32_t val = nd[i];
            nd[i] = (val >> 9) + carry;
            carry = (val & 0x1ff) * 1953125;
            if (i == (ndlo & 127)) break;
            i = (i - 1) & 127;
        }

        if (carry) nd[--ndlo & 127] = carry;
        k -= 9;
    }

    if (k) { // 0 =< k < 9 
        uint32_t mask = (1U << k) - 1;
        uint32_t mul = 1000000000 >> k;
        uint32_t i = ndhi & 127;
        uint32_t carry = 0;

        for (;;) {
            uint32_t val = nd[i];
            nd[i] = (val >> k) + carry;
            carry = (val & mask) * mul;
            if (i == (ndlo & 127)) break;
            i = (i - 1) & 127;
        }

        if (carry) nd[--ndlo & 127] = carry;
    }

    return ndlo;
}

// We can then go through the same process for multiplying, starting with 
// multiplication by 2:
// carry_in_max = val / 10^9 
//              = 4_294_967_295 / 10^9 = 4
// val_max = (10^9 - 1) * 2 + carry_in_max
//         = 2_000_000_002 < 4_294_967_295 (max of uint32_t)
int32_t nd_mul2(uint32_t* nd, int32_t ndhi) {
    uint32_t carry_in = 0;

    for (uint32_t i = 0; i <= (uint32_t)ndhi; i++) {
        uint32_t val = (nd[i] << 1) | carry_in;
        carry_in = val / 1000000000;
        // val = k * 10^9 + r (0 =< r < 10^9)
        // carry_in = k
        // nd[i] = val - carry_in * 10^9 = r
        nd[i] = val - carry_in * 1000000000;
    }

    if (carry_in) nd[++ndhi] = carry_in;

    return ndhi;
}

// By promoting val to 64-bits, this can be generalised to small k:
int32_t nd_mul2k_29(uint32_t* nd, int32_t ndhi, uint32_t k) {
    uint32_t carry_in = 0;

    for (uint32_t i = 0; i <= (uint32_t)ndhi; i++) {
        uint64_t val = ((uint64_t)nd[i] << k) | carry_in;
        carry_in = (uint32_t)(val / 1000000000);
        // val = k * 10^9 + r (0 =< r < 10^9)
        // carry_in = k
        // nd[i] = val - carry_in * 10^9 = r
        nd[i] = (uint32_t)(val - carry_in * 1000000000);
    }

    if (carry_in) nd[++ndhi] = carry_in;

    return ndhi;
}

// This time the constraint on k comes from wanting val to be no more
// than (10^9)^2 = 10^18, which limits k to 29. 
// carry_in_max = 10^9 - 1
// if k = 29
// val = (10^9 - 1) * 2^29 + carry_in_max
//     = 10^9 * 2^29 - 2^29 + 10^9 - 1
//     = 10^9 * 536_870_912 + 463,129,087
//     = 536_870_912_463_129_087 < 10^18
// if k = 30
// val = (10^9 - 1) * 2^30 + carry_in_max
//     = 10^9 * 2^30 - 2^30 + 10^9 - 1
//     = 1_073_741_824_000_000_000 - 73_741_825
//     = 1_073_741_823_926_258_174 > 10^18
// It also turns out to be useful to make carry_in a parameter, all of which 
// leads to the complete code for multiplying by 2^k:
int32_t nd_mul2k(uint32_t* nd, int32_t ndhi, uint32_t k, uint32_t carry_in) {

    while (k >= 29) {
        for (uint32_t i = 0; i <= (uint32_t)ndhi; i++) {
            uint64_t val = ((uint64_t)nd[i] << 29) | carry_in;
            carry_in = (uint32_t)(val / 1000000000);
            // val = k * 10^9 + r (0 =< r < 10^9)
            // carry_in = k
            // nd[i] = val - carry_in * 10^9 = r
            nd[i] = (uint32_t)(val - carry_in * 1000000000);
        }

        if (carry_in) {
            nd[++ndhi] = carry_in; 
            carry_in = 0;
        }
        k -= 29;
    }

    if (k) { // 0 =< k < 29 
        for (uint32_t i = 0; i <= (uint32_t)ndhi; i++) {
            uint64_t val = ((uint64_t)nd[i] << k) | carry_in;
            carry_in = (uint32_t)(val / 1000000000);
            // val = k * 10^9 + r (0 =< r < 10^9)
            // carry_in = k
            // nd[i] = val - carry_in * 10^9 = r
            nd[i] = (uint32_t)(val - carry_in * 1000000000);
        }

        if (carry_in) nd[++ndhi] = carry_in;
    }

    return ndhi;
}

// We can plug these routines into the decode function from earlier to create
// a print function:
void print(double n) {
    TValue t;
    t.n = n;
    if ((t.u32.hi << 1) >= 0xffe00000) {
        // binary form: 
        // t.u32.hi = Bx-111-1111-1111-xxxx-xxxx-xxxx-xxxx-xxxx
        // t.u32.hi << 1 = B1111-1111-111x-xxxx-xxxx-xxxx-xxxx-xxxx
        if (((t.u32.hi & 0x000fffff) | t.u32.lo) != 0) {
            printf("NaN\n");
        }
        else {
            printf("Infinity\n");
        }
    }
    else {
        char buf[1154];
        uint32_t nd[128];
        int32_t ndlo = 0;
        int32_t ndhi = 0;
        int32_t e = (t.u32.hi >> 20) & 0x7ff;
        nd[0] = t.u32.hi & 0xfffff;

        if (e == 0) {
            e++;
        }
        else {
            nd[0] |= 0x100000;
        }
        // max of nd[0] is 0x111111 = 2_097_151 < 10^9 ????

        e -= 1043;
        // if (t.u32.lo == 0)
        //     n = (m * 2^32 + t.u32.lo) * 2^(e - 1075) 
        //       = m * 2^32 * 2^(e - 1075) = m * 2^(e - 1043)
        // else 
        //     n = (m * 2^32 + t.u32.lo) * 2^(e - 1075)
        //       = (m * 2^3 + t.u32.lo / 2^29) * 2^29 * 2^(e - 1075) ????
        if (t.u32.lo) {
            e -= 32;
            // max of nd[0] is 0x111111000 +  = 2_097_151 < 10^9 ????
            nd[0] = (nd[0] << 3) | (t.u32.lo >> 29);
            ndhi = nd_mul2k(nd, ndhi, 29, t.u32.lo & 0x1fffffff);
        }

        if (e >= 0) {
            ndhi = nd_mul2k(nd, ndhi, (uint32_t)e, 0);
        }
        else {
            ndlo = nd_div2k(nd, ndlo, ndhi, (uint32_t)-e);
        }

        nd_print(buf, nd, ndlo, ndhi);
        printf("%s\n", buf);
    }
}