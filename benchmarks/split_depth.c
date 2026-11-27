/* split_depth: how deep pffft splits a transform (the number of split
 * levels its setup chooses, pffft_priv_impl.h) at given split thresholds,
 * for benchmarks/split_sim.py. A transform's speed depends on the threshold
 * only through this depth, so the time of a threshold that was not measured
 * is that of a measured one with the same depth.
 *
 *   split_depth f32|f64 T1,T2,... N ...
 *
 * prints "N real D1 D2 ..." and "N cplx D1 D2 ...", Di the depth at
 * threshold Ti (bytes, as split_max_real and split_max_complex), -1 where N
 * is not a valid length. The depth is the same on every SIMD backend (all
 * have 4-element vectors); without SIMD nothing is split. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int split_depth_f32(int n, int cplx, size_t split_max_real, size_t split_max_complex);
int split_depth_f64(int n, int cplx, size_t split_max_real, size_t split_max_complex);

int main(int argc, char **argv) {
  size_t thresholds[64];
  int nt = 0, f64, i, cplx, k;
  char *list, *tok;
  if (argc < 3 || (strcmp(argv[1], "f32") && strcmp(argv[1], "f64"))) {
    fprintf(stderr, "usage: split_depth f32|f64 T1,T2,... N ...\n");
    return 2;
  }
  f64 = !strcmp(argv[1], "f64");
  list = argv[2];
  for (tok = strtok(list, ","); tok && nt < 64; tok = strtok(NULL, ","))
    thresholds[nt++] = (size_t)strtoull(tok, NULL, 10);
  for (i = 3; i < argc; ++i) {
    int n = atoi(argv[i]);
    for (cplx = 0; cplx < 2; ++cplx) {
      printf("%d %s", n, cplx ? "cplx" : "real");
      for (k = 0; k < nt; ++k)
        printf(" %d", (f64 ? split_depth_f64 : split_depth_f32)(n, cplx, thresholds[k], thresholds[k]));
      printf("\n");
    }
  }
  return 0;
}
