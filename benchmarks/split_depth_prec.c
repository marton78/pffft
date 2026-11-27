/* One precision of split_depth (split_depth.c): the library's own sources
 * with the setup's fields in sight. Compiled once per precision, with
 * SPLIT_DEPTH_F64 for double. */
#ifdef SPLIT_DEPTH_F64
#  include "../src/pffft_double.c"
#  define SPLIT_DEPTH_FN split_depth_f64
#  define SPLIT_NEW pffftd_new_setup_tuned
#  define SPLIT_DEFAULT pffftd_default_tuning
#  define SPLIT_DESTROY pffftd_destroy_setup
#else
#  include "../src/pffft.c"
#  define SPLIT_DEPTH_FN split_depth_f32
#  define SPLIT_NEW pffft_new_setup_tuned
#  define SPLIT_DEFAULT pffft_default_tuning
#  define SPLIT_DESTROY pffft_destroy_setup
#endif

/* the split levels of a setup of n points, -1 if n is not valid */
int SPLIT_DEPTH_FN(int n, int cplx, size_t split_max_real, size_t split_max_complex) {
  pffft_tuning_t t;
  SETUP_STRUCT *s;
  int depth;
  SPLIT_DEFAULT(&t, sizeof t);
  t.split_max_real = split_max_real;
  t.split_max_complex = split_max_complex;
  s = SPLIT_NEW(n, cplx ? PFFFT_COMPLEX : PFFFT_REAL, &t);
  if (!s) return -1;
  depth = s->split;
  SPLIT_DESTROY(s);
  return depth;
}
