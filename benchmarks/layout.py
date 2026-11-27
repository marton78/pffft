#!/usr/bin/env python3
"""Where the code of a pffft library lies, and what it is.

    python3 benchmarks/layout.py LIB [LIB ...]
    python3 benchmarks/layout.py loops LIB [FUNCTION ...]
    python3 benchmarks/layout.py diff LIB_A LIB_B [-v FUNCTION ...]

The first form prints each library's functions with their address and its
offset modulo 64 (the alignment that decoders and op caches see).
compare_many.py --pads uses the functions here to show that each pad moved
the code as intended, and saves the addresses with its CSV.

loops lists, per function, where it starts and where its loops start (the
targets of backward branches inside it), with their offsets modulo 64 and
32 and their length in bytes.

diff compares the code of two libraries function by function, with the
addresses taken out (branch targets as FUNCTION+OFFSET, other addresses as
ADDR): which functions differ, by how many instructions. Two pads of one
revision must come out the same; for two revisions it shows which code
the change touched. -v prints the differences of the functions named.

The disassembly is llvm-objdump's (objdump on macOS is llvm-objdump; it reads
ELF, Mach-O and the DLLs that MSVC links, for x86 and Arm alike).

The addresses come from nm (llvm-nm or GNU nm: ELF and Mach-O), or, for a
library linked by MSVC, from the map file that layout_pad.cmake makes the
linker write (pffft.map in the build directory). Static functions are
included: they are most of pffft's kernels. A name that occurs more than once
(the same static function in the f32 and the f64 object) gets #2, #3, ... in
the order of the addresses.
"""
from __future__ import annotations

import bisect
import difflib
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path


class Error(Exception):
    pass


def text_symbols(lib: Path) -> dict:
    """{name: address} of the functions in lib, in address order."""
    lib = Path(lib)
    if lib.suffix.lower() == ".dll":
        found = map_symbols(find_map(lib))
    else:
        found = nm_symbols(lib)
    out, seen = {}, {}
    for addr, name in sorted(found):
        seen[name] = seen.get(name, 0) + 1
        out[name if seen[name] == 1 else f"{name}#{seen[name]}"] = addr
    return out


def nm_symbols(lib: Path) -> list:
    nm = next((t for t in ("llvm-nm", "nm") if shutil.which(t)), None)
    if nm is None:
        raise Error("no nm on the PATH")
    p = subprocess.run([nm, "--defined-only", str(lib)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode:
        raise Error(f"{nm} {lib}: {p.stderr.strip()}")
    out = []
    for line in p.stdout.splitlines():
        f = line.split()
        if len(f) == 3 and f[1] in "tT":
            name = f[2]
            if sys.platform == "darwin" or lib.suffix == ".dylib":
                name = name.removeprefix("_")  # Mach-O C names
            out.append((int(f[0], 16), name))
    return out


def find_map(lib: Path) -> Path:
    """LIB.map next to lib (a copy that compare_many.py saved), or pffft.map
    in the build directory of lib (the one with CMakeCache.txt)."""
    if lib.with_suffix(".map").is_file():
        return lib.with_suffix(".map")
    for d in [lib.parent, *lib.parents]:
        if (d / "pffft.map").is_file():
            return d / "pffft.map"
        if (d / "CMakeCache.txt").is_file():
            break
    raise Error(f"no pffft.map for {lib}: was it built with layout_pad.cmake?")


# " 0001:00000a40       passf4_ps                  0000000180001a40 f   pffft.obj"
MAP_LINE = re.compile(r"^\s*([0-9a-fA-F]{4}):[0-9a-fA-F]{8}\s+(\S+)\s+([0-9a-fA-F]{8,16})\s+(f\s+)?")


def map_symbols(path: Path) -> list:
    """Functions of an MSVC map file, publics and statics: the lines marked f
    (function); the pad, which is data, by its name."""
    out = []
    for line in path.read_text(errors="replace").splitlines():
        m = MAP_LINE.match(line)
        if m and (m.group(4) or m.group(2) == "pffft_layout_pad"):
            out.append((int(m.group(3), 16), m.group(2)))
    return out


def shift(base: dict, moved: dict) -> dict:
    """How far moved's functions lie from base's, for the names in both:
    {"shift": the most common difference, "same": how many functions moved by
    it, "of": of how many, "groups": {difference: functions}}. The objects
    of a library can move by different amounts (each by its pad, and by the
    pads of the objects before it), and the start-up code in front of the
    first object does not move."""
    common = [n for n in base if n in moved and not n.startswith("pffft_layout_pad")]
    if not common:
        return {"shift": None, "same": 0, "of": 0, "groups": {}}
    d = [moved[n] - base[n] for n in common]
    groups = {}
    for x in d:
        groups[x] = groups.get(x, 0) + 1
    mode = max(groups, key=groups.get)
    return {"shift": mode, "same": groups[mode], "of": len(d),
            "groups": dict(sorted(groups.items(), key=lambda kv: -kv[1]))}


def objdump() -> str:
    for tool in ("llvm-objdump", "objdump"):
        path = shutil.which(tool)
        if path:
            p = subprocess.run([path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if "LLVM" in p.stdout:
                return path
    raise Error("no llvm-objdump on the PATH (objdump on macOS is one)")


LINE = re.compile(r"^\s*([0-9a-fA-F]+):\s+(.*)$")
HEX = re.compile(r"0x([0-9a-fA-F]+)")
BRANCH = re.compile(r"^(j[a-z]*|loop[a-z]*|b(\.[a-z]+)?|b[a-z]{2}(\.w)?|cbn?z|tbn?z)$")


def disassemble(lib: Path) -> dict:
    """{function: [(address, instruction text), ...]} for the functions of
    text_symbols(lib), the comments of the disassembler taken out."""
    syms = text_symbols(lib)
    starts = sorted((a, n) for n, a in syms.items() if n != "pffft_layout_pad")
    addrs = [a for a, _ in starts]
    p = subprocess.run([objdump(), "-d", "--no-show-raw-insn", str(lib)], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, text=True)
    if p.returncode:
        raise Error(f"objdump {lib}: {p.stderr.strip()}")
    out, text_section = {}, False
    for line in p.stdout.splitlines():
        if line.startswith("Disassembly of section "):
            # only the code of the functions: not the stubs, the PLT, ...
            sect = line[len("Disassembly of section "):].rstrip(":").strip()
            text_section = sect in (".text", "__TEXT,__text")
            continue
        m = LINE.match(line)
        if not m or not text_section:
            continue
        addr = int(m.group(1), 16)
        k = bisect.bisect_right(addrs, addr) - 1
        if k < 0:
            continue
        text = re.split(r"\s+(?:;|//)\s|\s+#\s", m.group(2).strip())[0].strip()
        if not text or text.startswith("<"):
            continue
        out.setdefault(starts[k][1], []).append((addr, text))
    # Without symbols in the file (MSVC's DLLs) objdump decodes on from the
    # padding between functions, which can swallow a function's first bytes:
    # such a function is decoded again from its own start.
    ends = {n: (addrs[i + 1] if i + 1 < len(addrs) else None) for i, (_, n) in enumerate(starts)}
    for addr, name in starts:
        if name in out and out[name][0][0] != addr:
            out[name] = decode_range(lib, addr, ends[name]) or out[name]
    return out


def decode_range(lib: Path, start: int, stop) -> list:
    """[(address, text)] from start to stop (None: to the end of the section)."""
    cmd = [objdump(), "-d", "--no-show-raw-insn", f"--start-address={start:#x}", str(lib)]
    if stop is not None:
        cmd.insert(-1, f"--stop-address={stop:#x}")
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out = []
    for line in p.stdout.splitlines():
        m = LINE.match(line)
        if m:
            text = re.split(r"\s+(?:;|//)\s|\s+#\s", m.group(2).strip())[0].strip()
            if text and not text.startswith("<"):
                out.append((int(m.group(1), 16), text))
    return out


def normalize(insns, code: dict, immediates: bool = True) -> list:
    """The instructions without addresses: targets inside a function of code
    (disassemble()) as <FUNCTION+OFFSET>, other addresses (0x1000 and up)
    and displacements from the instruction pointer as ADDR. Immediates (#...
    on Arm, $... on x86) are kept, or with immediates=False replaced by IMM:
    on Arm the low bits of a data address are an immediate, which moves
    with the code."""
    starts = sorted((v[0][0], n) for n, v in code.items())
    addrs = [a for a, _ in starts]
    lo, hi = (addrs[0], max(v[-1][0] for v in code.values()) + 16) if addrs else (0, 0)

    def sub(m):
        v = int(m.group(2), 16)
        if lo <= v < hi:
            k = bisect.bisect_right(addrs, v) - 1
            return f"{m.group(1)}<{starts[k][1]}+{v - starts[k][0]:#x}>"
        return m.group(1) + ("ADDR" if v >= 0x1000 else m.group(0)[len(m.group(1)):])

    def outside(m):  # a target outside the functions: a stub, the PLT (the label of a stub can be wrong)
        v = int(m.group(1), 16)
        return m.group(0) if lo <= v < hi else "<EXT>"

    out = []
    for _, text in insns:
        text = re.sub(r"0x([0-9a-fA-F]+) <([^>]*)>", outside, text)
        # the disassembler's label after an address inside the functions goes: sub() names it
        text = " ".join(re.sub(r"(0x[0-9a-fA-F]+)\s*<[^>]*>", r"\1", text).split())
        if text.split()[0] == "adrp":  # a page: of data, or of anything
            text = re.sub(r",\s*\S+$", ", PAGE", text)
        text = re.sub(r"-?0x[0-9a-fA-F]+\(%rip\)", "ADDR(%rip)", text)
        text = re.sub(r"\[rip ?[+-] ?0x[0-9a-fA-F]+\]", "[rip+ADDR]", text)
        if not immediates:
            text = re.sub(r"[#$]-?(0x[0-9a-fA-F]+|\d+)", "IMM", text)
        out.append(re.sub(r"(^|[^#$\w])0x([0-9a-fA-F]+)", sub, text))
    return out


def loops(insns) -> list:
    """(head, end) of each backward branch inside the function: the loop
    runs from head to the branch at end."""
    if not insns:
        return []
    start, out = insns[0][0], []
    for addr, text in insns:
        mnemonic = text.split()[0].lower()
        if not BRANCH.match(mnemonic):
            continue
        for m in HEX.finditer(text):
            v = int(m.group(1), 16)
            if start <= v <= addr:
                out.append((v, addr))
                break
    return sorted(set(out))


def cmd_loops(lib: Path, names) -> None:
    code = disassemble(lib)
    for name, insns in code.items():
        if names and name not in names and name.split("#")[0] not in names:
            continue
        start = insns[0][0]
        print(f"{name}: {start:#x} (mod 64: {start % 64}), {insns[-1][0] - start + 4} bytes")
        for head, end in loops(insns):
            print(f"  loop +{head - start:#06x}: mod 64 {head % 64:2d}, mod 32 {head % 32:2d}, {end - head + 4:4d} bytes")


def cmd_diff(a: Path, b: Path, verbose) -> None:
    code_a, code_b = disassemble(a), disassemble(b)
    same = imm = 0
    for name in sorted(set(code_a) | set(code_b)):
        if name not in code_a or name not in code_b:
            print(f"only in {'B' if name in code_b else 'A'}: {name}")
            continue
        na, nb = normalize(code_a[name], code_a), normalize(code_b[name], code_b)
        if na == nb:
            same += 1
            continue
        if normalize(code_a[name], code_a, False) == normalize(code_b[name], code_b, False):
            imm += 1
            continue
        changed = sum(1 for op in difflib.SequenceMatcher(None, na, nb, autojunk=False).get_opcodes() if op[0] != "equal")
        print(f"differs: {name}: {len(na)} -> {len(nb)} instructions, {changed} changed runs")
        if name in verbose or name.split("#")[0] in verbose:
            sys.stdout.writelines(difflib.unified_diff([x + "\n" for x in na], [x + "\n" for x in nb], "A", "B", n=2))
    print(f"{same} functions the same, {imm} the same but for immediates (e.g. the low bits of data addresses on Arm)")


def main() -> None:
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        sys.exit(__doc__)
    if args[0] == "loops":
        cmd_loops(Path(args[1]), args[2:])
    elif args[0] == "diff":
        verbose = args[args.index("-v") + 1:] if "-v" in args else []
        cmd_diff(Path(args[1]), Path(args[2]), verbose)
    else:
        for lib in args:
            print(lib)
            for name, addr in text_symbols(Path(lib)).items():
                print(f"  {addr:#10x}  {addr % 64:2d}  {name}")


if __name__ == "__main__":
    try:
        main()
    except Error as e:
        sys.exit(f"layout.py: {e}")
