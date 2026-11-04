/* QEMU mps3-an547 (Cortex-M55): the vector table in ITCM at address 0,
   everything else in DDR, which leaves room for large C++ programs */
MEMORY { ITCM (rx) : ORIGIN = 0x00000000, LENGTH = 512K
         RAM (rwx) : ORIGIN = 0x60000000, LENGTH = 64M }
ENTRY(Reset_Handler)
SECTIONS {
  .vectors : { KEEP(*(.vectors)) } > ITCM
  .text : { *(.text*) KEEP(*(.init)) KEEP(*(.fini)) *(.rodata*)
            . = ALIGN(4); __preinit_array_start = .; KEEP(*(.preinit_array)) __preinit_array_end = .;
            __init_array_start = .; KEEP(*(SORT(.init_array.*))) KEEP(*(.init_array)) __init_array_end = .;
            __fini_array_start = .; KEEP(*(SORT(.fini_array.*))) KEEP(*(.fini_array)) __fini_array_end = .;
            KEEP(*(.eh_frame*)) } > RAM
  .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } > RAM
  .ARM.exidx : { __exidx_start = .; *(.ARM.exidx* .gnu.linkonce.armexidx.*) __exidx_end = .; } > RAM
  .data : { *(.data*) . = ALIGN(4); } > RAM
  .bss (NOLOAD) : { __bss_start__ = .; *(.bss*) *(COMMON) . = ALIGN(8); __bss_end__ = .; } > RAM
  end = .; _end = .; __end__ = .;
  __stack = ORIGIN(RAM) + LENGTH(RAM);
}
