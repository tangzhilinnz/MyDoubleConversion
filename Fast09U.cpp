#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>

#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

/*
 * Given an integer u from 0 to 9999, we want to perform 3 divisions
 * by constants 10, 100 and 1000 in parallel and calculate four digits
 * u - u/10*10, u/10 - u/100*10, etc. These digits can be shuffled,
 * converted to ascii and stored in memory as four consecutive bytes.
 *
 * One common approach to constant division is double-width multiplication
 * by a magic constant and shifting high-word to the right by a constant
 * number of bits.
 *
 * Double-width multiplication in xmm register can be done with pmuludq
 * but it operates on two 32-bit words while we need at least three
 * multiplications. For u that fits into 16-bit word, we can try pmaddwd
 * which multiplies eight signed 16-bit words, takes sums of pairs and
 * stores the results in four 32-bit words.
 *
 * The algorithm below uses these magic multiplications:
 *
 * u/10   : u * 26215 / 2^18,
 * u/100  : u * 10486 / 2^20,
 * u/1000 : u * 8389  / 2^23.
 *
 * The shifts are all different but it doesn't matter. Instead of
 * shifting to the right, low bits are masked and values are later
 * multiplied to scale the results by 256.
 */

void nasonov9(char* p, uint32_t u) {
    // __m128i _mm_set_epi16 (short e7, short e6, short e5, short e4, short e3,
    //                        short e2, short e1, short e0)
    // Set packed 16-bit integers in dst with the supplied values.
    // dst[15:0] := e0
    // dst[31:16] := e1
    // dst[47:32] := e2
    // dst[63:48] := e3
    // dst[79:64] := e4
    // dst[95:80] := e5
    // dst[111:96] := e6
    // dst[127:112] := e7

    // __m128i _mm_set1_epi16 (short a)
    // Broadcast 16-bit integer a to all all elements of dst. This intrinsic
    // may generate vpbroadcastw.
    // FOR j := 0 to 7
	// i := j*16
	// dst[i+15:i] := a[15:0]
    // ENDFOR

    // __m128i _mm_madd_epi16 (__m128i a, __m128i b)
    // Multiply packed signed 16-bit integers in a and b, producing intermediate
    // signed 32-bit integers. Horizontally add adjacent pairs of intermediate 
    // 32-bit integers, and pack the results in dst.
    // FOR j := 0 to 3
	// i := j*32
	// dst[i+31:i] := SignExtend32(a[i+31:i+16]*b[i+31:i+16]) + 
    //                SignExtend32(a[i+15:i]*b[i+15:i])
    // ENDFOR

    // __m128i _mm_and_si128 (__m128i a, __m128i b)
    // Compute the bitwise AND of 128 bits (representing integer data) in a and
    // b, and store the result in dst.
    // dst[127:0] := (a[127:0] AND b[127:0])

    // __m128i _mm_slli_si128 (__m128i a, int imm8)
    // Shift a left by imm8 bytes while shifting in zeros, and store the results
    // in dst.

    // __m128i _mm_or_si128 (__m128i a, __m128i b)
    // Compute the bitwise OR of 128 bits (representing integer data) in a and b,
    // and store the result in dst.
    // dst[127:0] := (a[127:0] OR b[127:0])

    // __m128i _mm_packs_epi32 (__m128i a, __m128i b)
    // Convert packed signed 32-bit integers from a and b to packed 16-bit integers
    // using signed saturation, and store the results in dst.
    // dst[15:0] := Saturate16(a[31:0])
    // dst[31:16] := Saturate16(a[63:32])
    // dst[47:32] := Saturate16(a[95:64])
    // dst[63:48] := Saturate16(a[127:96])
    // dst[79:64] := Saturate16(b[31:0])
    // dst[95:80] := Saturate16(b[63:32])
    // dst[111:96] := Saturate16(b[95:64])
    // dst[127:112] := Saturate16(b[127:96])

    // __m128i _mm_srli_epi16 (__m128i a, int imm8)
    // Shift packed 16-bit integers in a right by imm8 while shifting in zeros, 
    // and store the results in dst.

    // __m128i _mm_packs_epi16 (__m128i a, __m128i b)
    // Convert packed signed 16-bit integers from a and b 
    // to packed 8-bit integers using signed saturation, 
    // and store the results in dst.
    // dst[7:0] := Saturate8(a[15:0])
    // dst[15:8] := Saturate8(a[31:16])
    // dst[23:16] := Saturate8(a[47:32])
    // dst[31:24] := Saturate8(a[63:48])
    // dst[39:32] := Saturate8(a[79:64])
    // dst[47:40] := Saturate8(a[95:80])
    // dst[55:48] := Saturate8(a[111:96])
    // dst[63:56] := Saturate8(a[127:112])
    // dst[71:64] := Saturate8(b[15:0])
    // dst[79:72] := Saturate8(b[31:16])
    // dst[87:80] := Saturate8(b[47:32])
    // dst[95:88] := Saturate8(b[63:48])
    // dst[103:96] := Saturate8(b[79:64])
    // dst[111:104] := Saturate8(b[95:80])
    // dst[119:112] := Saturate8(b[111:96])
    // dst[127:120] := Saturate8(b[127:112])

    // void _mm_storel_epi64 (__m128i* mem_addr, __m128i a)
    // Store 64-bit integer from the first element of a into memory.
    // MEM[mem_addr+63:mem_addr] := a[63:0]

    uint32_t v = u / 10000;
    uint32_t w = v / 10000;
    u -= v * 10000;
    v -= w * 10000;

    const __m128i first_madd =
        _mm_set_epi16(-32768, -32768, 0, 26215, 0, 10486, 0, 8389);
    const __m128i mask =
        _mm_set_epi16(0xffff, 0, 0xfffc, 0, 0xfff0, 0, 0xff80, 0);
    const __m128i second_madd =
        _mm_set_epi16(-256, -640, 64, -160, 16, -20, 2, 0);

    // first_madd = [8389, 0, 10486, 0, 26215, 0, -32768, -32768]
    // mask = [0, 0xff80, 0, 0xfff0, 0, 0xfffc, 0, 0xffff]
    // second_madd = [0, 2, -20, 16, -160, 64, -640, -256]
    
    // short integer can perfectly contain v and u which are from 0 to 9999 
    // x = _mm_madd_epi16([v, v, v, v, v, v, v, v], first_madd)
    //   = [v * 8389, v * 10486, v * 26215, (v * -32768) + (v * -32768)]
    //   = [v * 8389, v * 10486, v * 26215, [0x0000, -v]] 
    // y = [u * 8389, u * 10486, u * 26215, [0x0000, -u]]
    __m128i x = _mm_madd_epi16(_mm_set1_epi16(v), first_madd);
    __m128i y = _mm_madd_epi16(_mm_set1_epi16(u), first_madd);

    // x = [0, v/1000*128, 0, v/100*16, 0, v/10*4, 0, -v]
    // y = [0, u/1000*128, 0, u/100*16, 0, u/10*4, 0, -u]
    x = _mm_and_si128(x, mask);
    y = _mm_and_si128(y, mask);

    // x = [0, v/1000*128, v/1000*128, v/100*16, v/100*16, v/10*4, v/10*4, -v]
	// y = [0, u/1000*128, u/1000*128, u/100*16, u/100*16, u/10*4, u/10*4, -u]
    x = _mm_or_si128(x, _mm_slli_si128(x, 2));
    y = _mm_or_si128(y, _mm_slli_si128(y, 2));

    // is multiplied to produce 4 scaled digits:
    // x = [(v/1000*128)*2, 
    //      (v/100*16)*16 - (v/1000*128)*2*10, 
    //      (v/10*4)*64 -(v/100*16)*16*10, 
    //      (-v)*-256 - (v/10*4)*10*64]
    // y = [(u/1000*128)*2, 
    //      (u/100*16)*16 - (u/1000*128)*2*10, 
    //      (u/10*4)*64 -(u/100*16)*16*10, 
    //      (-u)*-256 - (u/10*4)*10*64]
    x = _mm_madd_epi16(x, second_madd);
    y = _mm_madd_epi16(y, second_madd);

    // z = [v/1000, 
    //      v/100 -  v/1000 * 10, 
    //      v/10  -  v/100 * 10, 
    //      v     -  v/10 * 10,
    //      u/1000,
    //      u/100 - u/1000 * 10, 
    //      u/10  - u/100 * 10, 
    //      u     - u/10 * 10]
    __m128i z = _mm_srli_epi16(_mm_packs_epi32(x, y), 8);

    // z = [v/1000, 
    //      v/100 -  v/1000 * 10, 
    //      v/10  -  v/100 * 10, 
    //      v     -  v/10 * 10,
    //      u/1000,
    //      u/100 - u/1000 * 10, 
    //      u/10  - u/100 * 10, 
    //      u     - u/10 * 10,
    //
    //      v/1000, 
    //      v/100 -  v/1000 * 10, 
    //      v/10  -  v/100 * 10, 
    //      v     -  v/10 * 10,
    //      u/1000,
    //      u/100 - u/1000 * 10, 
    //      u/10  - u/100 * 10, 
    //      u     - u/10 * 10]
    z = _mm_packs_epi16(z, z);
    p[0] = '0' | w;
    _mm_storel_epi64((__m128i*)(p + 1), _mm_or_si128(z, _mm_set1_epi32(0x30303030)));
}

void divmod9(char* p, uint32_t n) {
    for (uint32_t i = 9; i--; n /= 10) {
        p[i] = '0' + (n % 10);
    }
}


int main() {

    auto start = system_clock::now();

    char buf[100] = { 0 };

    for (int i = 0; i < 100000000; i++) {
        nasonov9(buf, i);
		//divmod9(buf, i);
    }

    std::cout << "Hello World!\n" << buf << std::endl;

    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);

    std::cout << "cost: "
        << double(duration.count()) * microseconds::period::num / microseconds::period::den << "seconds" << std::endl;
}