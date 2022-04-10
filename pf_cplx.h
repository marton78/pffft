/*
This software is part of pffft/pfdsp, a set of simple DSP routines.

Copyright (c) 2020  Hayati Ayguen <h_ayguen@web.de>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif


/*
   _____                      _
  / ____|                    | |
 | |     ___  _ __ ___  _ __ | | _____  __
 | |    / _ \| '_ ` _ \| '_ \| |/ _ \ \/ /
 | |___| (_) | | | | | | |_) | |  __/>  <
  \_____\___/|_| |_| |_| .__/|_|\___/_/\_\
                       | |
                       |_|
*/

typedef struct complexf_s {
    float i;
    float q;

#ifdef __cplusplus
    inline complexf_s & operator=(float r) {
        this->i = r;
        this->q = 0.0F;
        return *this;
    }
#endif

} complexf;


typedef struct complexd_s {
    double i;
    double q;

#ifdef __cplusplus
    inline complexd_s & operator=(double r) {
        this->i = r;
        this->q = 0.0;
        return *this;
    }
#endif

} complexd;



#ifdef __cplusplus
}

/* add/define some helper functions */

static inline float real(const float self) {
    return self;
}

static inline float real(const complexf &self) {
    return self.i;
}

static inline float imag(const complexf &self) {
    return self.q;
}

static inline void operator+=(complexf &self, float r) {
    self.i += r;
}

static inline void operator+=(complexf &self, const complexf sum) {
    self.i += sum.i;
    self.q += sum.q;
}

static inline void operator*=(complexf &self, float r) {
    self.i *= r;
    self.q *= r;
}



static inline double real(const double self) {
    return self;
}

static inline double real(const complexd &self) {
    return self.i;
}

static inline double imag(const complexd &self) {
    return self.q;
}

static inline void operator+=(complexd &self, double r) {
    self.i += r;
}

static inline void operator+=(complexd &self, const complexd sum) {
    self.i += sum.i;
    self.q += sum.q;
}

static inline void operator+=(complexd &self, const complexf sum) {
    self.i += sum.i;
    self.q += sum.q;
}

static inline void operator*=(complexd &self, double r) {
    self.i *= r;
    self.q *= r;
}

#endif
