#ifndef RecoTracker_MkFitCore_src_Matriplex_MatriplexSym_h
#define RecoTracker_MkFitCore_src_Matriplex_MatriplexSym_h

#include "MatriplexCommon.h"
#include "Matriplex.h"

//==============================================================================
// MatriplexSym
//==============================================================================

namespace Matriplex {

  const idx_t gSymOffsets[7][36] = {
      {},
      {},
      {0, 1, 1, 2},
      {0, 1, 3, 1, 2, 4, 3, 4, 5},  // 3
      {},
      {0, 1, 3, 6, 10, 1, 2, 4, 7, 11, 3, 4, 5, 8, 12, 6, 7, 8, 9, 13, 10, 11, 12, 13, 14},
      {0, 1, 3, 6, 10, 15, 1,  2,  4,  7,  11, 16, 3,  4,  5,  8,  12, 17,
       6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20}};

  //------------------------------------------------------------------------------

  template <typename T, idx_t D, idx_t N>
  class __attribute__((aligned(MPLEX_ALIGN))) MatriplexSym {
  public:
    typedef T value_type;

    /// no. of matrix rows
    static constexpr int kRows = D;
    /// no. of matrix columns
    static constexpr int kCols = D;
    /// no of elements: lower triangle
    static constexpr int kSize = (D + 1) * D / 2;
    /// size of the whole matriplex
    static constexpr int kTotSize = N * kSize;

    T fArray[kTotSize];

    MatriplexSym() {}
    MatriplexSym(T v) { setVal(v); }

    idx_t plexSize() const { return N; }

    void setVal(T v) {
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] = v;
      }
    }

    void add(const MatriplexSym& v) {
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] += v.fArray[i];
      }
    }

    void scale(T scale) {
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] *= scale;
      }
    }

    MatriplexSym& negate() {
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] = -fArray[i];
      }
      return *this;
    }

    template <typename TT>
    MatriplexSym& negate_if_ltz(const MatriplexSym<TT, D, N>& sign) {
      for (idx_t i = 0; i < kTotSize; ++i) {
        if (sign.fArray[i] < 0)
          fArray[i] = -fArray[i];
      }
      return *this;
    }

    T operator[](idx_t xx) const { return fArray[xx]; }
    T& operator[](idx_t xx) { return fArray[xx]; }

    const idx_t* offsets() const { return gSymOffsets[D]; }
    idx_t off(idx_t i) const { return gSymOffsets[D][i]; }

    const T& constAt(idx_t n, idx_t i, idx_t j) const { return fArray[off(i * D + j) * N + n]; }

    T& At(idx_t n, idx_t i, idx_t j) { return fArray[off(i * D + j) * N + n]; }

    T& operator()(idx_t n, idx_t i, idx_t j) { return At(n, i, j); }
    const T& operator()(idx_t n, idx_t i, idx_t j) const { return constAt(n, i, j); }

    // reduction/assignment helpers

    using QReduced = Matriplex<T, 1, 1, N>;

    QReduced ReduceFixedIJ(idx_t i, idx_t j) const {
      QReduced t;
      for (idx_t n = 0; n < N; ++n) {
        t[n] = constAt(n, i, j);
      }
      return t;
    }
    QReduced rij(idx_t i, idx_t j) const { return ReduceFixedIJ(i, j); }
    QReduced operator()(idx_t i, idx_t j) const { return ReduceFixedIJ(i, j); }

    struct QAssigner {
      MatriplexSym& m_matriplex;
      const int m_i, m_j;

      QAssigner(MatriplexSym& m, int i, int j) : m_matriplex(m), m_i(i), m_j(j) {}
      MatriplexSym& operator=(const QReduced& qvec) {
        for (idx_t n = 0; n < N; ++n) {
          m_matriplex(n, m_i, m_j) = qvec[n];
        }
        return m_matriplex;
      }
      MatriplexSym& operator=(T qscalar) {
        for (idx_t n = 0; n < N; ++n) {
          m_matriplex(n, m_i, m_j) = qscalar;
        }
        return m_matriplex;
      }
    };

    QAssigner AssignFixedIJ(idx_t i, idx_t j) { return QAssigner(*this, i, j); }
    QAssigner aij(idx_t i, idx_t j) { return AssignFixedIJ(i, j); }

    MatriplexSym& operator=(const MatriplexSym& m) {
      memcpy(fArray, m.fArray, sizeof(T) * kTotSize);
      return *this;
    }

    MatriplexSym(const MatriplexSym& m) = default;

    MatriplexSym& operator=(T t) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = t;
      return *this;
    }

    MatriplexSym& operator+=(T t) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] += t;
      return *this;
    }

    MatriplexSym& operator-=(T t) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] -= t;
      return *this;
    }

    MatriplexSym& operator*=(T t) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] *= t;
      return *this;
    }

    MatriplexSym& operator/=(T t) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] /= t;
      return *this;
    }

    MatriplexSym& operator+=(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] += a.fArray[i];
      return *this;
    }

    MatriplexSym& operator-=(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] -= a.fArray[i];
      return *this;
    }

    MatriplexSym& operator*=(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] *= a.fArray[i];
      return *this;
    }

    MatriplexSym& operator/=(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] /= a.fArray[i];
      return *this;
    }

    MatriplexSym operator-() {
      MatriplexSym t;
      for (idx_t i = 0; i < kTotSize; ++i)
        t.fArray[i] = -fArray[i];
      return t;
    }

    MatriplexSym& abs(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::abs(a.fArray[i]);
      return *this;
    }
    MatriplexSym& abs() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::abs(fArray[i]);
      return *this;
    }

    MatriplexSym& sqr(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = a.fArray[i] * a.fArray[i];
      return *this;
    }
    MatriplexSym& sqr() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = fArray[i] * fArray[i];
      return *this;
    }

    MatriplexSym& sqrt(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::sqrt(a.fArray[i]);
      return *this;
    }
    MatriplexSym& sqrt() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::sqrt(fArray[i]);
      return *this;
    }

    MatriplexSym& hypot(const MatriplexSym& a, const MatriplexSym& b) {
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] = a.fArray[i] * a.fArray[i] + b.fArray[i] * b.fArray[i];
      }
      return sqrt();
    }

    MatriplexSym& sin(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::sin(a.fArray[i]);
      return *this;
    }
    MatriplexSym& sin() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::sin(fArray[i]);
      return *this;
    }

    MatriplexSym& cos(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::cos(a.fArray[i]);
      return *this;
    }
    MatriplexSym& cos() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::cos(fArray[i]);
      return *this;
    }

    MatriplexSym& tan(const MatriplexSym& a) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::tan(a.fArray[i]);
      return *this;
    }
    MatriplexSym& tan() {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::tan(fArray[i]);
      return *this;
    }

    MatriplexSym& atan2(const MatriplexSym& y, const MatriplexSym& x) {
      for (idx_t i = 0; i < kTotSize; ++i)
        fArray[i] = std::atan2(y.fArray[i], x.fArray[i]);
      return *this;
    }

    void sincos(MatriplexSym& s, MatriplexSym& c) const {
      for (idx_t i = 0; i < kTotSize; ++i) {
        s.fArray[i] = std::sin(fArray[i]);
        c.fArray[i] = std::cos(fArray[i]);
      }
    }

    void sincos4(MatriplexSym& s, MatriplexSym& c) const {
      for (idx_t i = 0; i < kTotSize; ++i)
        internal::sincos4(fArray[i], s.fArray[i], c.fArray[i]);
    }

#ifdef MPLEX_VDT

#define ASS fArray[i] =
#define ARR fArray[i]
#define A_ARR a.fArray[i]

#ifdef MPLEX_VDT_USE_STD
#define VDT_INVOKE(_ass_, _func_, ...) \
  for (idx_t i = 0; i < kTotSize; ++i) \
    _ass_ std::_func_(__VA_ARGS__);
#else
#define VDT_INVOKE(_ass_, _func_, ...)          \
  for (idx_t i = 0; i < kTotSize; ++i)          \
    if constexpr (std::is_same<T, float>())     \
      _ass_ vdt::fast_##_func_##f(__VA_ARGS__); \
    else                                        \
      _ass_ vdt::fast_##_func_(__VA_ARGS__);
#endif

    MatriplexSym& fast_isqrt(const MatriplexSym& a) {
      VDT_INVOKE(ASS, isqrt, A_ARR);
      return *this;
    }
    MatriplexSym& fast_isqrt() {
      VDT_INVOKE(ASS, isqrt, ARR);
      return *this;
    }

    MatriplexSym& fast_sin(const MatriplexSym& a) {
      VDT_INVOKE(ASS, sin, A_ARR);
      return *this;
    }
    MatriplexSym& fast_sin() {
      VDT_INVOKE(ASS, sin, ARR);
      return *this;
    }

    MatriplexSym& fast_cos(const MatriplexSym& a) {
      VDT_INVOKE(ASS, cos, A_ARR);
      return *this;
    }
    MatriplexSym& fast_cos() {
      VDT_INVOKE(ASS, cos, ARR);
      return *this;
    }

    void fast_sincos(MatriplexSym& s, MatriplexSym& c) const { VDT_INVOKE(, sincos, ARR, s.fArray[i], c.fArray[i]); }

    MatriplexSym& fast_tan(const MatriplexSym& a) {
      VDT_INVOKE(ASS, tan, A_ARR);
      return *this;
    }
    MatriplexSym& fast_tan() {
      VDT_INVOKE(ASS, tan, ARR);
      return *this;
    }

    MatriplexSym& fast_atan2(const MatriplexSym& y, const MatriplexSym& x) {
      VDT_INVOKE(ASS, atan2, y.fArray[i], x.fArray[i]);
      return *this;
    }

#undef VDT_INVOKE

#undef ASS
#undef ARR
#undef A_ARR
#endif

    void copySlot(idx_t n, const MatriplexSym& m) {
      for (idx_t i = n; i < kTotSize; i += N) {
        fArray[i] = m.fArray[i];
      }
    }

    void copyIn(idx_t n, const T* arr) {
      for (idx_t i = n; i < kTotSize; i += N) {
        fArray[i] = *(arr++);
      }
    }

    void copyIn(idx_t n, const MatriplexSym& m, idx_t in) {
      for (idx_t i = n; i < kTotSize; i += N, in += N) {
        fArray[i] = m[in];
      }
    }

    void copy(idx_t n, idx_t in) {
      for (idx_t i = n; i < kTotSize; i += N, in += N) {
        fArray[i] = fArray[in];
      }
    }

#if defined(AVX512_INTRINSICS)

    template <typename U>
    void slurpIn(const T* arr, __m512i& vi, const U&, const int N_proc = N) {
      //_mm512_prefetch_i32gather_ps(vi, arr, 1, _MM_HINT_T0);

      const __m512 src = {0};
      const __mmask16 k = N_proc == N ? -1 : (1 << N_proc) - 1;

      for (int i = 0; i < kSize; ++i, ++arr) {
        //_mm512_prefetch_i32gather_ps(vi, arr+2, 1, _MM_HINT_NTA);

        __m512 reg = _mm512_mask_i32gather_ps(src, k, vi, arr, sizeof(U));
        _mm512_mask_store_ps(&fArray[i * N], k, reg);
      }
    }

    // Experimental methods, slurpIn() seems to be at least as fast.
    // See comments in mkFit/MkFitter.cc MkFitter::addBestHit().

    void ChewIn(const char* arr, int off, int vi[N], const char* tmp, __m512i& ui) {
      // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.

      for (int i = 0; i < N; ++i) {
        __m512 reg = _mm512_load_ps(arr + vi[i]);
        _mm512_store_ps((void*)(tmp + 64 * i), reg);
      }

      for (int i = 0; i < kSize; ++i) {
        __m512 reg = _mm512_i32gather_ps(ui, tmp + off + i * sizeof(T), 1);
        _mm512_store_ps(&fArray[i * N], reg);
      }
    }

    void Contaginate(const char* arr, int vi[N], const char* tmp) {
      // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.

      for (int i = 0; i < N; ++i) {
        __m512 reg = _mm512_load_ps(arr + vi[i]);
        _mm512_store_ps((void*)(tmp + 64 * i), reg);
      }
    }

    void Plexify(const char* tmp, __m512i& ui) {
      for (int i = 0; i < kSize; ++i) {
        __m512 reg = _mm512_i32gather_ps(ui, tmp + i * sizeof(T), 1);
        _mm512_store_ps(&fArray[i * N], reg);
      }
    }

#elif defined(AVX2_INTRINSICS)

    template <typename U>
    void slurpIn(const T* arr, __m256i& vi, const U&, const int N_proc = N) {
      const __m256 src = {0};

      __m256i k = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
      __m256i k_sel = _mm256_set1_epi32(N_proc);
      __m256i k_master = _mm256_cmpgt_epi32(k_sel, k);

      k = k_master;
      for (int i = 0; i < kSize; ++i, ++arr) {
        __m256 reg = _mm256_mask_i32gather_ps(src, arr, vi, (__m256)k, sizeof(U));
        // Restore mask (docs say gather clears it but it doesn't seem to).
        k = k_master;
        _mm256_maskstore_ps(&fArray[i * N], k, reg);
      }
    }

#else

    void slurpIn(const T* arr, int vi[N], const int N_proc = N) {
      // Separate N_proc == N case (gains about 7% in fit test).
      if (N_proc == N) {
        for (int i = 0; i < kSize; ++i) {
          for (int j = 0; j < N; ++j) {
            fArray[i * N + j] = *(arr + i + vi[j]);
          }
        }
      } else {
        for (int i = 0; i < kSize; ++i) {
          for (int j = 0; j < N_proc; ++j) {
            fArray[i * N + j] = *(arr + i + vi[j]);
          }
        }
      }
    }

#endif

    void copyOut(idx_t n, T* arr) const {
      for (idx_t i = n; i < kTotSize; i += N) {
        *(arr++) = fArray[i];
      }
    }

    void setDiagonal3x3(idx_t n, T d) {
      T* p = fArray + n;

      p[0 * N] = d;
      p[1 * N] = 0;
      p[2 * N] = d;
      p[3 * N] = 0;
      p[4 * N] = 0;
      p[5 * N] = d;
    }

    MatriplexSym& subtract(const MatriplexSym& a, const MatriplexSym& b) {
      // Does *this = a - b;

#pragma omp simd
      for (idx_t i = 0; i < kTotSize; ++i) {
        fArray[i] = a.fArray[i] - b.fArray[i];
      }

      return *this;
    }

    // ==================================================================
    // Operations specific to Kalman fit in 6 parameter space
    // ==================================================================

    void addNoiseIntoUpperLeft3x3(T noise) {
      T* p = fArray;
      ASSUME_ALIGNED(p, 64);

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
        p[0 * N + n] += noise;
        p[2 * N + n] += noise;
        p[5 * N + n] += noise;
      }
    }

    void invertUpperLeft3x3() {
      typedef T TT;

      T* a = fArray;
      ASSUME_ALIGNED(a, 64);

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
        const TT c00 = a[2 * N + n] * a[5 * N + n] - a[4 * N + n] * a[4 * N + n];
        const TT c01 = a[4 * N + n] * a[3 * N + n] - a[1 * N + n] * a[5 * N + n];
        const TT c02 = a[1 * N + n] * a[4 * N + n] - a[2 * N + n] * a[3 * N + n];
        const TT c11 = a[5 * N + n] * a[0 * N + n] - a[3 * N + n] * a[3 * N + n];
        const TT c12 = a[3 * N + n] * a[1 * N + n] - a[4 * N + n] * a[0 * N + n];
        const TT c22 = a[0 * N + n] * a[2 * N + n] - a[1 * N + n] * a[1 * N + n];

        // Force determinant calculation in double precision.
        const double det = (double)a[0 * N + n] * c00 + (double)a[1 * N + n] * c01 + (double)a[3 * N + n] * c02;
        const TT s = TT(1) / det;

        a[0 * N + n] = s * c00;
        a[1 * N + n] = s * c01;
        a[2 * N + n] = s * c11;
        a[3 * N + n] = s * c02;
        a[4 * N + n] = s * c12;
        a[5 * N + n] = s * c22;
      }
    }

  };

  template <typename T, idx_t D, idx_t N>
  using MPlexSym = MatriplexSym<T, D, N>;

  //==============================================================================
  // Operators
  //==============================================================================

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator-(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t = a;
    t.negate();
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> negate(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t = a;
    t.negate();
    return t;
  }

  template <typename T, typename TT, idx_t D, idx_t N>
  MPlexSym<T, D, N> negate_if_ltz(const MPlexSym<T, D, N>& a, const MPlexSym<TT, D, N>& sign) {
    MPlexSym<T, D, N> t = a;
    t.negate_if_ltz(sign);
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator+(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t += b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator-(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t -= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator*(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t *= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator/(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t /= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator+(const MPlexSym<T, D, N>& a, T b) {
    MPlexSym<T, D, N> t = a;
    t += b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator-(const MPlexSym<T, D, N>& a, T b) {
    MPlexSym<T, D, N> t = a;
    t -= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator*(const MPlexSym<T, D, N>& a, T b) {
    MPlexSym<T, D, N> t = a;
    t *= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator/(const MPlexSym<T, D, N>& a, T b) {
    MPlexSym<T, D, N> t = a;
    t /= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator+(T a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t += b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator-(T a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t -= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator*(T a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t *= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> operator/(T a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t = a;
    t /= b;
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> abs(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.abs(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> sqr(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.sqr(a);
  }

  //---------------------------------------------------------
  // transcendentals, std version

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> sqrt(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.sqrt(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> hypot(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t;
    return t.hypot(a, b);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> sin(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.sin(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> cos(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.cos(a);
  }

  template <typename T, idx_t D, idx_t N>
  void sincos(const MPlexSym<T, D, N>& a, MPlexSym<T, D, N>& s, MPlexSym<T, D, N>& c) {
    a.sincos(s, c);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> tan(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.tan(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> atan2(const MPlexSym<T, D, N>& y, const MPlexSym<T, D, N>& x) {
    MPlexSym<T, D, N> t;
    return t.atan2(y, x);
  }

  //---------------------------------------------------------
  // transcendentals, vdt version

#ifdef MPLEX_VDT

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> fast_isqrt(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.fast_isqrt(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> fast_sin(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.fast_sin(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> fast_cos(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.fast_cos(a);
  }

  template <typename T, idx_t D, idx_t N>
  void fast_sincos(const MPlexSym<T, D, N>& a, MPlexSym<T, D, N>& s, MPlexSym<T, D, N>& c) {
    a.fast_sincos(s, c);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> fast_tan(const MPlexSym<T, D, N>& a) {
    MPlexSym<T, D, N> t;
    return t.fast_tan(a);
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> fast_atan2(const MPlexSym<T, D, N>& y, const MPlexSym<T, D, N>& x) {
    MPlexSym<T, D, N> t;
    return t.fast_atan2(y, x);
  }

#endif

  template <typename T, idx_t D, idx_t N>
  void sincos4(const MPlexSym<T, D, N>& a, MPlexSym<T, D, N>& s, MPlexSym<T, D, N>& c) {
    a.sincos4(s, c);
  }

  //---------------------------------------------------------

  template <typename T, idx_t D, idx_t N>
  void min_max(const MPlexSym<T, D, N>& a,
               const MPlexSym<T, D, N>& b,
               MPlexSym<T, D, N>& min,
               MPlexSym<T, D, N>& max) {
    for (idx_t i = 0; i < a.kTotSize; ++i) {
      min.fArray[i] = std::min(a.fArray[i], b.fArray[i]);
      max.fArray[i] = std::max(a.fArray[i], b.fArray[i]);
    }
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> min(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t;
    for (idx_t i = 0; i < a.kTotSize; ++i) {
      t.fArray[i] = std::min(a.fArray[i], b.fArray[i]);
    }
    return t;
  }

  template <typename T, idx_t D, idx_t N>
  MPlexSym<T, D, N> max(const MPlexSym<T, D, N>& a, const MPlexSym<T, D, N>& b) {
    MPlexSym<T, D, N> t;
    for (idx_t i = 0; i < a.kTotSize; ++i) {
      t.fArray[i] = std::max(a.fArray[i], b.fArray[i]);
    }
    return t;
  }

  //==============================================================================
  // Multiplications
  //==============================================================================

  template <typename T, idx_t D, idx_t N>
  struct SymMultiplyCls {
    static void multiply(const MPlexSym<T, D, N>& A, const MPlexSym<T, D, N>& B, MPlex<T, D, D, N>& C) {
      throw std::runtime_error("general symmetric multiplication not supported");
    }
  };

  template <typename T, idx_t N>
  struct SymMultiplyCls<T, 3, N> {
    static void multiply(const MPlexSym<T, 3, N>& A, const MPlexSym<T, 3, N>& B, MPlex<T, 3, 3, N>& C) {
      const T* a = A.fArray;
      ASSUME_ALIGNED(a, 64);
      const T* b = B.fArray;
      ASSUME_ALIGNED(b, 64);
      T* c = C.fArray;
      ASSUME_ALIGNED(c, 64);

#ifdef MPLEX_INTRINSICS

      for (idx_t n = 0; n < N; n += 64 / sizeof(T)) {
#include "intr_sym_3x3.ah"
      }

#else

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
#include "std_sym_3x3.ah"
      }

#endif
    }
  };

  template <typename T, idx_t N>
  struct SymMultiplyCls<T, 6, N> {
    static void multiply(const MPlexSym<float, 6, N>& A, const MPlexSym<float, 6, N>& B, MPlex<float, 6, 6, N>& C) {
      const T* a = A.fArray;
      ASSUME_ALIGNED(a, 64);
      const T* b = B.fArray;
      ASSUME_ALIGNED(b, 64);
      T* c = C.fArray;
      ASSUME_ALIGNED(c, 64);

#ifdef MPLEX_INTRINSICS

      for (idx_t n = 0; n < N; n += 64 / sizeof(T)) {
#include "intr_sym_6x6.ah"
      }

#else

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
#include "std_sym_6x6.ah"
      }

#endif
    }
  };

  template <typename T, idx_t D, idx_t N>
  void multiply(const MPlexSym<T, D, N>& A, const MPlexSym<T, D, N>& B, MPlex<T, D, D, N>& C) {
    SymMultiplyCls<T, D, N>::multiply(A, B, C);
  }

  //==============================================================================
  // Cramer inversion
  //==============================================================================

  template <typename T, idx_t D, idx_t N>
  struct CramerInverterSym {
    static void invert(MPlexSym<T, D, N>& A, double* determ = nullptr) {
      throw std::runtime_error("general cramer inversion not supported");
    }
  };

  template <typename T, idx_t N>
  struct CramerInverterSym<T, 2, N> {
    static void invert(MPlexSym<T, 2, N>& A, double* determ = nullptr) {
      typedef T TT;

      T* a = A.fArray;
      ASSUME_ALIGNED(a, 64);

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
        // Force determinant calculation in double precision.
        const double det = (double)a[0 * N + n] * a[2 * N + n] - (double)a[1 * N + n] * a[1 * N + n];
        if (determ)
          determ[n] = det;

        const TT s = TT(1) / det;
        const TT tmp = s * a[2 * N + n];
        a[1 * N + n] *= -s;
        a[2 * N + n] = s * a[0 * N + n];
        a[0 * N + n] = tmp;
      }
    }
  };

  template <typename T, idx_t N>
  struct CramerInverterSym<T, 3, N> {
    static void invert(MPlexSym<T, 3, N>& A, double* determ = nullptr) {
      typedef T TT;

      T* a = A.fArray;
      ASSUME_ALIGNED(a, 64);

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
        const TT c00 = a[2 * N + n] * a[5 * N + n] - a[4 * N + n] * a[4 * N + n];
        const TT c01 = a[4 * N + n] * a[3 * N + n] - a[1 * N + n] * a[5 * N + n];
        const TT c02 = a[1 * N + n] * a[4 * N + n] - a[2 * N + n] * a[3 * N + n];
        const TT c11 = a[5 * N + n] * a[0 * N + n] - a[3 * N + n] * a[3 * N + n];
        const TT c12 = a[3 * N + n] * a[1 * N + n] - a[4 * N + n] * a[0 * N + n];
        const TT c22 = a[0 * N + n] * a[2 * N + n] - a[1 * N + n] * a[1 * N + n];

        // Force determinant calculation in double precision.
        const double det = (double)a[0 * N + n] * c00 + (double)a[1 * N + n] * c01 + (double)a[3 * N + n] * c02;
        if (determ)
          determ[n] = det;

        const TT s = TT(1) / det;
        a[0 * N + n] = s * c00;
        a[1 * N + n] = s * c01;
        a[2 * N + n] = s * c11;
        a[3 * N + n] = s * c02;
        a[4 * N + n] = s * c12;
        a[5 * N + n] = s * c22;
      }
    }
  };

  template <typename T, idx_t D, idx_t N>
  void invertCramerSym(MPlexSym<T, D, N>& A, double* determ = nullptr) {
    CramerInverterSym<T, D, N>::invert(A, determ);
  }

  //==============================================================================
  // Cholesky inversion
  //==============================================================================

  template <typename T, idx_t D, idx_t N>
  struct CholeskyInverterSym {
    static void invert(MPlexSym<T, D, N>& A) { throw std::runtime_error("general cholesky inversion not supported"); }
  };

  template <typename T, idx_t N>
  struct CholeskyInverterSym<T, 3, N> {
    static void invert(MPlexSym<T, 3, N>& A) {
      typedef T TT;

      T* a = A.fArray;

#pragma omp simd
      for (idx_t n = 0; n < N; ++n) {
        TT l0 = std::sqrt(T(1) / a[0 * N + n]);
        TT l1 = a[1 * N + n] * l0;
        TT l2 = a[2 * N + n] - l1 * l1;
        l2 = std::sqrt(T(1) / l2);
        TT l3 = a[3 * N + n] * l0;
        TT l4 = (a[4 * N + n] - l1 * l3) * l2;
        TT l5 = a[5 * N + n] - (l3 * l3 + l4 * l4);
        l5 = std::sqrt(T(1) / l5);

        // decomposition done

        l3 = (l1 * l4 * l2 - l3) * l0 * l5;
        l1 = -l1 * l0 * l2;
        l4 = -l4 * l2 * l5;

        a[0 * N + n] = l3 * l3 + l1 * l1 + l0 * l0;
        a[1 * N + n] = l3 * l4 + l1 * l2;
        a[2 * N + n] = l4 * l4 + l2 * l2;
        a[3 * N + n] = l3 * l5;
        a[4 * N + n] = l4 * l5;
        a[5 * N + n] = l5 * l5;

        // m(2,x) are all zero if anything went wrong at l5.
        // all zero, if anything went wrong already for l0 or l2.
      }
    }
  };

  template <typename T, idx_t D, idx_t N>
  void invertCholeskySym(MPlexSym<T, D, N>& A) {
    CholeskyInverterSym<T, D, N>::invert(A);
  }

}  // end namespace Matriplex

#endif
