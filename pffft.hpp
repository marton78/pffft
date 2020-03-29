/* Copyright (c) 2020  Dario Mambro ( dario.mambro@gmail.com )
   Copyright (c) 2020  Hayati Ayguen ( h_ayguen@web.de )

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of PFFFT, nor the names of its
   sponsors or contributors may be used to endorse or promote products
   derived from this Software without specific prior written permission.

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
   SOFTWARE.
*/

#pragma once

#include <complex>

namespace {
#include "pffft.h"
#include "pffft_double.h"
}

namespace pffft {

// enum { PFFFT_REAL, PFFFT_COMPLEX }
typedef pffft_transform_t TransformType;

// enum { PFFFT_FORWARD, PFFFT_BACKWARD }
typedef pffft_direction_t TransformDirection;

namespace {
template<typename T>
class Setup;
}

// T can be float, double, std::complex<float> or std::complex<double>
template<typename T>
class Fft
{
public:
  typedef typename Setup<T>::Scalar Scalar;
  typedef T TType;

  Fft(int length, int stackThresholdLen = 4096);

  ~Fft();

  void prepareLength(int newLength);

  int getLength() const { return length; }

  std::complex<Scalar>* forward(const T* input, std::complex<Scalar>* spectrum);

  T* inverse(const std::complex<Scalar>* spectrum, T* output);

  Scalar* forwardInternalLayout(const T* input,
                                Scalar* spectrum_internal_layout);

  T* inverseInternalLayout(const Scalar* spectrum_internal_layout, T* output);

  void reorderSpectrum(const Scalar* input,
                       std::complex<Scalar>* output,
                       TransformDirection direction);

  Scalar* convolveAccumulate(const Scalar* dft_a,
                             const Scalar* dft_b,
                             Scalar* dft_ab,
                             const Scalar scaling);

  Scalar* convolve(const Scalar* dft_a,
                   const Scalar* dft_b,
                   Scalar* dft_ab,
                   const Scalar scaling);


  int getSpectrumSize() const { return isComplexTransform() ? length : ( length / 2 ); }

  int getInternalLayoutSize() const { return isComplexTransform() ? ( 2 * length ) : length; }


  T* allocateOrigin() const { return alignedAllocType(length); }

  std::complex<Scalar>* allocateSpectrum() const { return alignedAllocComplex( getSpectrumSize() ); }

  Scalar* allocateInternalLayout() const { return alignedAllocScalar( getInternalLayoutSize() ); }


  static bool isComplexTransform() { return sizeof(T) != sizeof(Scalar); }

  static bool isFloatScalar() { return sizeof(Scalar) == sizeof(float); }

  static bool isDoubleScalar() { return sizeof(Scalar) == sizeof(double); }

  static void alignedFree(void* ptr);

  static T * alignedAllocType(int length);

  static Scalar* alignedAllocScalar(int length);

  static std::complex<Scalar>* alignedAllocComplex(int length);

private:
  Setup<T> setup;
  Scalar* work;
  int length;
  int stackThresholdLen;
};

// helpers

void
alignedFree(void* ptr)
{
  pffft_aligned_free(ptr);
}

/* simple helper to get minimum possible fft length */
int
minFFtsize(const TransformType transform)
{
  return pffft_min_fft_size(transform);
}

/* simple helper to determine next power of 2
   - without inexact/rounding floating point operations
*/
int
nextPowerOfTwo(const int N)
{
  return pffft_next_power_of_two(N);
}

/* simple helper to determine if power of 2 - returns bool */
bool
isPowerOfTwo(const int N)
{
  return pffft_is_power_of_two(N);
}

// implementation

namespace {

template<typename T>
struct Setup
{};

template<>
class Setup<float>
{
  PFFFT_Setup* self;

public:
  typedef float Scalar;
  typedef float TType;

  Setup()
    : self(NULL)
  {}

  void prepareLength(int length)
  {
    if (self) {
      pffft_destroy_setup(self);
    }
    self = pffft_new_setup(length, PFFFT_REAL);
  }

  ~Setup() { pffft_destroy_setup(self); }

  void transform_ordered(const Scalar* input,
                         Scalar* output,
                         Scalar* work,
                         pffft_direction_t direction)
  {
    pffft_transform_ordered(self, input, output, work, direction);
  }

  void transform(const Scalar* input,
                 Scalar* output,
                 Scalar* work,
                 pffft_direction_t direction)
  {
    pffft_transform(self, input, output, work, direction);
  }

  void reorder(const Scalar* input, Scalar* output, pffft_direction_t direction)
  {
    pffft_zreorder(self, input, output, direction);
  }

  void convolveAccumulate(const Scalar* dft_a,
                          const Scalar* dft_b,
                          Scalar* dft_ab,
                          const Scalar scaling)
  {
    pffft_zconvolve_accumulate(self, dft_a, dft_b, dft_ab, scaling);
  }

  void convolve(const Scalar* dft_a,
                const Scalar* dft_b,
                Scalar* dft_ab,
                const Scalar scaling)
  {
    pffft_zconvolve_no_accu(self, dft_a, dft_b, dft_ab, scaling);
  }

  static TType* allocateType(int length)
  {
    const int bytes = sizeof(TType) * length;
    return static_cast<TType*>(pffft_aligned_malloc(bytes));
  }

  static Scalar* allocate(int length)
  {
    const int bytes = sizeof(Scalar) * length;
    return static_cast<Scalar*>(pffft_aligned_malloc(bytes));
  }
};

template<>
class Setup< std::complex<float> >
{
  PFFFT_Setup* self;

public:
  typedef float Scalar;
  typedef std::complex<float> TType;

  Setup()
    : self(NULL)
  {}

  ~Setup() { pffft_destroy_setup(self); }

  void prepareLength(int length)
  {
    if (self) {
      pffft_destroy_setup(self);
    }
    self = pffft_new_setup(length, PFFFT_COMPLEX);
  }

  void transform_ordered(const Scalar* input,
                         Scalar* output,
                         Scalar* work,
                         pffft_direction_t direction)
  {
    pffft_transform_ordered(self, input, output, work, direction);
  }

  void transform(const Scalar* input,
                 Scalar* output,
                 Scalar* work,
                 pffft_direction_t direction)
  {
    pffft_transform(self, input, output, work, direction);
  }

  void reorder(const Scalar* input, Scalar* output, pffft_direction_t direction)
  {
    pffft_zreorder(self, input, output, direction);
  }

  void convolve(const Scalar* dft_a,
                const Scalar* dft_b,
                Scalar* dft_ab,
                const Scalar scaling)
  {
    pffft_zconvolve_no_accu(self, dft_a, dft_b, dft_ab, scaling);
  }

  static TType* allocateType(int length)
  {
    const int bytes = sizeof(TType) * length;
    return static_cast<TType*>(pffft_aligned_malloc(bytes));
  }

  static Scalar* allocate(const int length)
  {
    const int bytes = sizeof(Scalar) * length;
    return static_cast<Scalar*>(pffft_aligned_malloc(bytes));
  }
};

template<>
class Setup<double>
{
  PFFFTD_Setup* self;

public:
  typedef double Scalar;
  typedef double TType;

  Setup()
    : self(NULL)
  {}

  ~Setup() { pffftd_destroy_setup(self); }

  void prepareLength(int length)
  {
    if (self) {
      pffftd_destroy_setup(self);
      self = NULL;
    }
    if (length > 0) {
      self = pffftd_new_setup(length, PFFFT_REAL);
    }
  }

  void transform_ordered(const Scalar* input,
                         Scalar* output,
                         Scalar* work,
                         pffft_direction_t direction)
  {
    pffftd_transform_ordered(self, input, output, work, direction);
  }

  void transform(const Scalar* input,
                 Scalar* output,
                 Scalar* work,
                 pffft_direction_t direction)
  {
    pffftd_transform(self, input, output, work, direction);
  }

  void reorder(const Scalar* input, Scalar* output, pffft_direction_t direction)
  {
    pffftd_zreorder(self, input, output, direction);
  }

  void convolveAccumulate(const Scalar* dft_a,
                          const Scalar* dft_b,
                          Scalar* dft_ab,
                          const Scalar scaling)
  {
    pffftd_zconvolve_accumulate(self, dft_a, dft_b, dft_ab, scaling);
  }

  void convolve(const Scalar* dft_a,
                const Scalar* dft_b,
                Scalar* dft_ab,
                const Scalar scaling)
  {
    pffftd_zconvolve_no_accu(self, dft_a, dft_b, dft_ab, scaling);
  }

  static TType* allocateType(int length)
  {
    const int bytes = sizeof(TType) * length;
    return static_cast<TType*>(pffft_aligned_malloc(bytes));
  }

  static Scalar* allocate(int length)
  {
    const int bytes = sizeof(Scalar) * length;
    return static_cast<Scalar*>(pffftd_aligned_malloc(bytes));
  }
};

template<>
class Setup< std::complex<double> >
{
  PFFFTD_Setup* self;

public:
  typedef double Scalar;
  typedef std::complex<double> TType;

  Setup()
    : self(NULL)
  {}

  ~Setup() { pffftd_destroy_setup(self); }

  void prepareLength(int length)
  {
    if (self) {
      pffftd_destroy_setup(self);
    }
    self = pffftd_new_setup(length, PFFFT_COMPLEX);
  }

  void transform_ordered(const Scalar* input,
                         Scalar* output,
                         Scalar* work,
                         pffft_direction_t direction)
  {
    pffftd_transform_ordered(self, input, output, work, direction);
  }

  void transform(const Scalar* input,
                 Scalar* output,
                 Scalar* work,
                 pffft_direction_t direction)
  {
    pffftd_transform(self, input, output, work, direction);
  }

  void reorder(const Scalar* input, Scalar* output, pffft_direction_t direction)
  {
    pffftd_zreorder(self, input, output, direction);
  }

  void convolveAccumulate(const Scalar* dft_a,
                          const Scalar* dft_b,
                          Scalar* dft_ab,
                          const Scalar scaling)
  {
    pffftd_zconvolve_accumulate(self, dft_a, dft_b, dft_ab, scaling);
  }

  void convolve(const Scalar* dft_a,
                const Scalar* dft_b,
                Scalar* dft_ab,
                const Scalar scaling)
  {
    pffftd_zconvolve_no_accu(self, dft_a, dft_b, dft_ab, scaling);
  }

  static TType* allocateType(int length)
  {
    const int bytes = sizeof(TType) * length;
    return static_cast<TType*>(pffft_aligned_malloc(bytes));
  }

  static Scalar* allocate(int length)
  {
    const int bytes = sizeof(Scalar) * length;
    return static_cast<Scalar*>(pffftd_aligned_malloc(bytes));
  }
};

}

template<typename T>
inline Fft<T>::Fft(int length, int stackThresholdLen)
  : length(0)
  , stackThresholdLen(stackThresholdLen)
  , work(NULL)
{
  prepareLength(length);
}

template<typename T>
inline Fft<T>::~Fft()
{
  pffft_aligned_free(work);
}

template<typename T>
inline void
Fft<T>::prepareLength(int newLength)
{
  const bool wasOnHeap = ( work != NULL );

  const bool useHeap = newLength > stackThresholdLen;

  if (useHeap == wasOnHeap && newLength == length) {
    return;
  }

  length = newLength;

  setup.prepareLength(length);

  if (work) {
    pffft_aligned_free(work);
    work = NULL;
  }

  if (useHeap) {
    int const bytesToAllocate = length * sizeof(T);
    work = static_cast<Scalar*>(pffft_aligned_malloc(bytesToAllocate));
  }
}

template<typename T>
inline std::complex<typename Fft<T>::Scalar>*
Fft<T>::forward(const T* input, std::complex<Scalar>* spectrum)
{
  setup.transform_ordered(reinterpret_cast<const Scalar*>(input),
                          reinterpret_cast<Scalar*>(spectrum),
                          work,
                          PFFFT_FORWARD);
  return spectrum;
}

template<typename T>
inline T*
Fft<T>::inverse(std::complex<Scalar> const* spectrum, T* output)
{
  setup.transform_ordered(reinterpret_cast<const Scalar*>(spectrum),
                          reinterpret_cast<Scalar*>(output),
                          work,
                          PFFFT_BACKWARD);
  return output;
}

template<typename T>
inline typename pffft::Fft<T>::Scalar*
Fft<T>::forwardInternalLayout(const T* input, Scalar* spectrum_internal_layout)
{
  setup.transform(reinterpret_cast<const Scalar*>(input),
                  spectrum_internal_layout,
                  work,
                  PFFFT_FORWARD);
  return spectrum_internal_layout;
}

template<typename T>
inline T*
Fft<T>::inverseInternalLayout(const Scalar* spectrum_internal_layout, T* output)
{
  setup.transform(spectrum_internal_layout,
                  reinterpret_cast<Scalar*>(output),
                  work,
                  PFFFT_BACKWARD);
  return output;
}

template<typename T>
inline void
Fft<T>::reorderSpectrum(const Scalar* input,
                        std::complex<Scalar>* output,
                        TransformDirection direction)
{
  setup.reorder(input, reinterpret_cast<Scalar*>(output), direction);
}

template<typename T>
inline typename pffft::Fft<T>::Scalar*
Fft<T>::convolveAccumulate(const Scalar* dft_a,
                           const Scalar* dft_b,
                           Scalar* dft_ab,
                           const Scalar scaling)
{
  setup.convolveAccumulate(dft_a, dft_b, dft_ab, scaling);
  return dft_ab;
}

template<typename T>
inline typename pffft::Fft<T>::Scalar*
Fft<T>::convolve(const Scalar* dft_a,
                 const Scalar* dft_b,
                 Scalar* dft_ab,
                 const Scalar scaling)
{
  setup.convolve(dft_a, dft_b, dft_ab, scaling);
  return dft_ab;
}

template<typename T>
inline void
Fft<T>::alignedFree(void* ptr)
{
  pffft_aligned_free(ptr);
}


template<typename T>
inline T*
pffft::Fft<T>::alignedAllocType(int length)
{
  return Setup<T>::allocateType(length);
}

template<typename T>
inline typename pffft::Fft<T>::Scalar*
pffft::Fft<T>::alignedAllocScalar(int length)
{
  return reinterpret_cast< Scalar* >( Setup<T>::allocate(length) );
}

template<typename T>
inline std::complex<typename pffft::Fft<T>::Scalar>*
Fft<T>::alignedAllocComplex(int length)
{
  return reinterpret_cast< std::complex<Scalar>* >( Setup<T>::allocate(2 * length) );
}

} // namespace pffft

