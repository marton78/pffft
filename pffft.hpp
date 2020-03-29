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
#include <vector>
#include <limits>

namespace {
#include "pffft.h"
#include "pffft_double.h"
}

namespace pffft {

// enum { PFFFT_REAL, PFFFT_COMPLEX }
typedef pffft_transform_t TransformType;

template<typename T>
class PFAlloc;

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
  typedef T value_type;

  typedef std::vector< T, PFAlloc<T> > ValueVector;
  typedef std::vector< std::complex<Scalar>, PFAlloc< std::complex<Scalar> > > ComplexVector;
  typedef std::vector< Scalar, PFAlloc<Scalar> > InternalLayoutVector;

  // static retrospection functions

  static bool isComplexTransform() { return sizeof(T) != sizeof(Scalar); }

  static bool isFloatScalar() { return sizeof(Scalar) == sizeof(float); }

  static bool isDoubleScalar() { return sizeof(Scalar) == sizeof(double); }

  //////////////////

  Fft( int length, int stackThresholdLen = 4096 );

  ~Fft();

  // set transformation length,
  //   which is identical to forward()'s input vector's size,
  //   and also equals inverse()'s output vector size
  void prepareLength(int newLength);

  int getLength() const { return length; }

  int getSpectrumSize() const { return isComplexTransform() ? length : ( length / 2 ); }

  int getInternalLayoutSize() const { return isComplexTransform() ? ( 2 * length ) : length; }


  ////////////////////////////////////////////
  ////
  //// API 1, with std::vector<> based containers,
  ////   which free the allocated memory themselves (RAII).
  ////
  //// uses an Allocator for the alignment of SIMD data.
  ////
  ////////////////////////////////////////////

  ValueVector valueVector() const { return ValueVector(length); }

  ComplexVector spectrumVector() const { return ComplexVector( getSpectrumSize() ); }

  InternalLayoutVector internalLayoutVector() const { return InternalLayoutVector( getInternalLayoutSize() ); }

  ////////////////////////////////////////////
  // although using Vectors for output ..
  // they need to have resize() applied before!

  // core API, having the spectrum in canonical order

  ComplexVector& forward(const ValueVector& input, ComplexVector& spectrum)
    { forward( input.data(), spectrum.data() ); return spectrum; }

  ValueVector& inverse(const ComplexVector& spectrum, ValueVector& output)
    { inverse( spectrum.data(), output.data() ); return output; }

  // provide additional functions with spectrum in some internal Layout.
  // these are faster, cause the implementation omits the reordering.
  // these are useful in special applications, like fast convolution,
  // where inverse() is following anyway ..

  InternalLayoutVector& forwardToInternalLayout(
          const ValueVector& input,
          InternalLayoutVector& spectrum_internal_layout
        )
    { forwardToInternalLayout( input.data(), spectrum_internal_layout.data() ); return spectrum_internal_layout; }

  ValueVector& inverseFromInternalLayout(
          const InternalLayoutVector& spectrum_internal_layout,
          ValueVector& output
        )
    { inverseFromInternalLayout( spectrum_internal_layout.data(), output.data() ); return output; }

  void reorderSpectrum(
          const InternalLayoutVector& input,
          ComplexVector &output
        )
    { reorderSpectrum( input.data(), output.data() ); }

  InternalLayoutVector& convolveAccumulate(
          const InternalLayoutVector& spectrum_internal_a,
          const InternalLayoutVector& spectrum_internal_b,
          InternalLayoutVector& spectrum_internal_ab,
          const Scalar scaling
        )
    { convolveAccumulate( spectrum_internal_a.data(), spectrum_internal_b.data(), spectrum_internal_ab.data(), scaling ); return spectrum_internal_ab; }

  InternalLayoutVector& convolve(
          const InternalLayoutVector& spectrum_internal_a,
          const InternalLayoutVector& spectrum_internal_b,
          InternalLayoutVector& spectrum_internal_ab,
          const Scalar scaling
        )
    { convolve( spectrum_internal_a.data(), spectrum_internal_b.data(), spectrum_internal_ab.data(), scaling ); return spectrum_internal_ab; }




  ////////////////////////////////////////////
  ////
  //// API 2, dealing with raw pointers,
  //// which need to be deallocated using alignedFree()
  ////
  //// the special allocation is required cause SIMD
  //// implementations require aligned memory
  ////
  ////////////////////////////////////////////

  static void alignedFree(void* ptr);

  static T * alignedAllocType(int length);

  static Scalar* alignedAllocScalar(int length);

  static std::complex<Scalar>* alignedAllocComplex(int length);

  // core API, having the spectrum in canonical order

  std::complex<Scalar>* forward(const T* input, std::complex<Scalar>* spectrum);

  T* inverse(const std::complex<Scalar>* spectrum, T* output);


  // provide additional functions with spectrum in some internal Layout.
  // these are faster, cause the implementation omits the reordering.
  // these are useful in special applications, like fast convolution,
  // where inverse() is following anyway ..

  Scalar* forwardToInternalLayout(const T* input,
                                Scalar* spectrum_internal_layout);

  T* inverseFromInternalLayout(const Scalar* spectrum_internal_layout, T* output);

  void reorderSpectrum(const Scalar* input, std::complex<Scalar>* output );

  Scalar* convolveAccumulate(const Scalar* spectrum_internal_a,
                             const Scalar* spectrum_internal_b,
                             Scalar* spectrum_internal_ab,
                             const Scalar scaling);

  Scalar* convolve(const Scalar* spectrum_internal_a,
                   const Scalar* spectrum_internal_b,
                   Scalar* spectrum_internal_ab,
                   const Scalar scaling);


private:
  Setup<T> setup;
  Scalar* work;
  int length;
  int stackThresholdLen;
};


////////////////////////////////////////////
////
//// PUBLIC HELPER FUNCTIONS
////
////////////////////////////////////////////

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

void
alignedFree(void* ptr)
{
  pffft_aligned_free(ptr);
}


////////////////////////////////////////////////////////////////////

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
  typedef float value_type;

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

  static value_type* allocateType(int length)
  {
    const int bytes = sizeof(value_type) * length;
    return static_cast<value_type*>(pffft_aligned_malloc(bytes));
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
  typedef std::complex<float> value_type;

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

  static value_type* allocateType(int length)
  {
    const int bytes = sizeof(value_type) * length;
    return static_cast<value_type*>(pffft_aligned_malloc(bytes));
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
  typedef double value_type;

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

  static value_type* allocateType(int length)
  {
    const int bytes = sizeof(value_type) * length;
    return static_cast<value_type*>(pffft_aligned_malloc(bytes));
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
  typedef std::complex<double> value_type;

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

  static value_type* allocateType(int length)
  {
    const int bytes = sizeof(value_type) * length;
    return static_cast<value_type*>(pffft_aligned_malloc(bytes));
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
Fft<T>::forwardToInternalLayout(const T* input, Scalar* spectrum_internal_layout)
{
  setup.transform(reinterpret_cast<const Scalar*>(input),
                  spectrum_internal_layout,
                  work,
                  PFFFT_FORWARD);
  return spectrum_internal_layout;
}

template<typename T>
inline T*
Fft<T>::inverseFromInternalLayout(const Scalar* spectrum_internal_layout, T* output)
{
  setup.transform(spectrum_internal_layout,
                  reinterpret_cast<Scalar*>(output),
                  work,
                  PFFFT_BACKWARD);
  return output;
}

template<typename T>
inline void
Fft<T>::reorderSpectrum( const Scalar* input, std::complex<Scalar>* output )
{
  setup.reorder(input, reinterpret_cast<Scalar*>(output), PFFFT_FORWARD);
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



////////////////////////////////////////////////////////////////////

// Allocator - for std::vector<>:
// origin: http://www.josuttis.com/cppcode/allocator.html
// http://www.josuttis.com/cppcode/myalloc.hpp
//
// minor renaming and utilizing of pffft (de)allocation functions
// are applied to Jossutis' allocator

/* The following code example is taken from the book
 * "The C++ Standard Library - A Tutorial and Reference"
 * by Nicolai M. Josuttis, Addison-Wesley, 1999
 *
 * (C) Copyright Nicolai M. Josuttis 1999.
 * Permission to copy, use, modify, sell and distribute this software
 * is granted provided this copyright notice appears in all copies.
 * This software is provided "as is" without express or implied
 * warranty, and with no claim as to its suitability for any purpose.
 */

template <class T>
class PFAlloc {
  public:
    // type definitions
    typedef T        value_type;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef T&       reference;
    typedef const T& const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    // rebind allocator to type U
    template <class U>
    struct rebind {
        typedef PFAlloc<U> other;
    };

    // return address of values
    pointer address (reference value) const {
        return &value;
    }
    const_pointer address (const_reference value) const {
        return &value;
    }

    /* constructors and destructor
     * - nothing to do because the allocator has no state
     */
    PFAlloc() throw() {
    }
    PFAlloc(const PFAlloc&) throw() {
    }
    template <class U>
      PFAlloc (const PFAlloc<U>&) throw() {
    }
    ~PFAlloc() throw() {
    }

    // return maximum number of elements that can be allocated
    size_type max_size () const throw() {
        return std::numeric_limits<std::size_t>::max() / sizeof(T);
    }

    // allocate but don't initialize num elements of type T
    pointer allocate (size_type num, const void* = 0) {
        pointer ret = (pointer)( pffft_aligned_malloc(num*sizeof(T)) );
        return ret;
    }

    // initialize elements of allocated storage p with value value
    void construct (pointer p, const T& value) {
        // initialize memory with placement new
        new((void*)p)T(value);
    }

    // destroy elements of initialized storage p
    void destroy (pointer p) {
        // destroy objects by calling their destructor
        p->~T();
    }

    // deallocate storage p of deleted elements
    void deallocate (pointer p, size_type num) {
        // deallocate memory with pffft
        pffft_aligned_free((void*)p);
    }
};

// return that all specializations of this allocator are interchangeable
template <class T1, class T2>
bool operator== (const PFAlloc<T1>&,
                 const PFAlloc<T2>&) throw() {
    return true;
}
template <class T1, class T2>
bool operator!= (const PFAlloc<T1>&,
                 const PFAlloc<T2>&) throw() {
    return false;
}


} // namespace pffft

