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
    } else {
      printf("Infinity\n");
    }
  } else {
    int32_t e = (t.u32.hi >> 20) & 0x7ff;
    uint64_t m = t.u32.hi & 0xfffff;
    if (e == 0) {
      e++;
    } else {
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
  uint32_t i = ndhi & 127, carry = 0;

  for (;;) {
    uint32_t val = nd[i];

    // (val & 1) is the first bit of the binary form of piece_i before 
	// dividing by 2, 
	// (nd[i] * (10^9)^i) / 2 
	//     = (val - (val & 1) + (val & 1)) / 2 * (10^9)^i
	//     = (val >> 1) * (10^9)^i + (val & 1) / 2 * (10^9)^i
	//     = (val >> 1) * (10^9)^i + (val & 1) * 500000000 * (10^9)^(i-1)
    nd[i] = (val >> 1) + carry;
    carry = (val & 1) * 500000000;
    if (i == (ndlo & 127)) break;
    i = (i - 1) & 127;
  }

  if (carry) nd[--ndlo & 127] = carry;
  return ndlo;
}

// We can generalise this to dividing by 2^k:
int32_t nd_div2k(uint32_t* nd, int32_t ndlo, int32_t ndhi, uint32_t k) {
  uint32_t mask = (1U << k) - 1, mul = 1000000000 >> k;
  uint32_t i = ndhi & 127, carry = 0;

  for (;;) {
    uint32_t val = nd[i];

	// (val & mask) is the first k bits of the binary form of piece_i before 
	// dividing by 2^k, 
	// (nd[i] * (10^9)^i) / 2^k 
	//     = (val - (val & mask) + (val & mask)) / 2^k * (10^9)^i
	//     = (val >> k) * (10^9)^i + (val & mask) / 2^k * (10^9)^i
	//     = (val >> k) * (10^9)^i + (val & mask) * ((10^9) >> k) * (10^9)^(i-1)
    nd[i] = (val >> k) + carry;
    carry = (val & mask) * mul;
    if (i == (ndlo & 127)) break;
    i = (i - 1) & 127;
  }

  if (carry) nd[--ndlo & 127] = carry;
  return ndlo;
}