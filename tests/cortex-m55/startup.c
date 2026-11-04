/* Minimal Cortex-M startup for running newlib/rdimon programs under QEMU. */
#include <errno.h>
#include <stddef.h>

extern void _start(void);           /* newlib's rdimon crt0 */
extern unsigned long __stack;       /* from the linker script */
extern char end;                    /* end of .bss, from the linker script */

static void sh_write0(const char *s) {
  register int r0 __asm__("r0") = 0x04; register const char *r1 __asm__("r1") = s;
  __asm__ volatile ("bkpt 0xab" : "+r"(r0) : "r"(r1) : "memory");
}
static void sh_exit_failure(void) {
  register int r0 __asm__("r0") = 0x18; register int r1 __asm__("r1") = 0x20023; /* RunTimeErrorUnknown */
  __asm__ volatile ("bkpt 0xab" : "+r"(r0) : "r"(r1) : "memory");
  for (;;) {}
}
void Fault_Handler(void) { sh_write0("*** CPU fault (HardFault/MemManage/BusFault/UsageFault)\n"); sh_exit_failure(); }
void Reset_Handler(void) {
  /* enable CP10/CP11 (FPU, and MVE on Armv8.1-M) before any FP instruction */
  *(volatile unsigned long *)0xE000ED88 |= (0xFUL << 20);
  __asm__ volatile ("dsb\n isb");
  _start();
}
__attribute__((section(".vectors"), used))
void (* const vectors[16])(void) = {
  (void (*)(void))&__stack, Reset_Handler,
  Fault_Handler, Fault_Handler, Fault_Handler, Fault_Handler, Fault_Handler,
};

/* Heap from the linker script, not from SYS_HEAPINFO, which QEMU guesses
   wrongly on this board: from the end of .bss to 256 KiB below the stack. */
void *_sbrk(ptrdiff_t incr) {
  static char *brk = &end;
  char *limit = (char *)&__stack - 256 * 1024;
  if (brk + incr > limit) { errno = ENOMEM; return (void *)-1; }
  char *prev = brk; brk += incr; return prev;
}
