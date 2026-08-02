# Full-cumulant 1-2-2 derivation — run state

Working state for the staged derivation (`tools/mkno_stage.wl`). Rebuildable:
delete `results/stage_*.m` and re-run each stage.

Goal: full-size 1-2-2 triple commutator with cumulants RETAINED (no
`NoCumulants`), then report term count, how many carry lambda, and the DELTA
against the NoCumulants-reduced version. Commit results and tell John.

Invocation: `MKNO_STAGE=<name> MKNO_MEMCAP=<GB> MKNO_TIMECAP=<s> tools/wl -file tools/mkno_stage.wl`

| stage | status | terms | with lambda | notes |
|---|---|---|---|---|
| r1v1  | built | 30 | 15 | <1 s, peak 71 MB |
| v1r1  | built | 30 | 15 | <1 s, peak 71 MB |
| rv1v2 | running (job 17ef0024-b677-462f-9ee3-c0755e5bb3c1) | — | — | 6 GB / 2400 s caps |
| v1rv2 | not started | | | |
| v2v1r | not started | | | |
| v2rv1 | not started | | | |
| combine | not started | | | emits mkno_full_cumulant.m + _nocum.m + census |

## History / hazards

- The unstaged run (`tools/mkno_full_cumulant.wl`, job 07941007) grew ~2 GB per
  30 s and was stopped at 3.0 GB before the OOM killer got it. Exit 143 is that
  deliberate stop, not a crash.
- **Never SIGKILL a Wolfram kernel.** It does not release its license lease, and
  the original entitlement allowed exactly one kernel, which locked the machine
  out with "No valid password found". Replacement entitlement has a 2-kernel
  limit and a 5-minute polling interval. Stop work with the in-script
  Memory/TimeConstrained guards instead.
- `wolframscript` launches the kernel under its own wrapper, so stopping the
  background job does NOT stop the kernel; it must be signalled directly.
- If `rv1v2` does not fit in 6 GB: combine identical monomials in the 30-term
  intermediates before multiplying. That is exact bookkeeping, not truncation,
  so it does not violate John's "take the size" instruction.

## Machine limits (measured 2026-08-02, not assumed)

`/proc/meminfo` reports **32 GB and 32 CPUs — that is the HOST**. /proc is not
namespaced in this container, so those numbers are misleading.

The real limit is the cgroup:

| knob | value |
|---|---|
| `/sys/fs/cgroup/memory.max` | **8 GB** — hard cap, process is OOM-killed above it |
| `/sys/fs/cgroup/memory.swap.max` | 8 GB allowed |
| host `SwapFree` | ~7 GB |
| **effective ceiling** | **~15 GB** (8 RAM + ~7 swap, the swap part slow) |

`/sys/fs/cgroup` is mounted **read-only** and `mount -o remount,rw` is denied
(no CAP_SYS_ADMIN in this namespace), so the 8 GB cap **cannot be raised from
inside the sandbox**. It is a platform allocation, not a configurable.

Consequence: `MKNO_MEMCAP` default raised 6 -> 12 GB in `tools/mkno_stage.wl`.
A stage larger than that must be reduced algebraically (combine like monomials
in the intermediates before multiplying -- exact bookkeeping, not truncation),
not given more memory.
