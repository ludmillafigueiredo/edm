slurmstepd: *** JOB 48585 ON uranus2 CANCELLED AT 2018-05-04T08:12:56 ***

signal (15): Beendet
while loading /home/luf74xx/Dokumente/model/main.jl, in expression starting on line 213
macro expansion at ./multidimensional.jl:527 [inlined]
macro expansion at ./cartesian.jl:64 [inlined]
macro expansion at ./multidimensional.jl:525 [inlined]
_unsafe_getindex! at ./multidimensional.jl:519 [inlined]
macro expansion at ./multidimensional.jl:513 [inlined]
_unsafe_getindex at ./multidimensional.jl:506
macro expansion at ./multidimensional.jl:495 [inlined]
_getindex at ./multidimensional.jl:491
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
jl_apply at /buildworker/worker/package_linux64/build/src/julia.h:1424 [inlined]
jl_f__apply at /buildworker/worker/package_linux64/build/src/builtins.c:426
projvegmass! at /home/luf74xx/Dokumente/model/Organisms.jl:128
unknown function (ip: 0x7f0495c270cd)
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
simulate at /home/luf74xx/Dokumente/model/main.jl:175
unknown function (ip: 0x7f049a9d59cf)
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
do_call at /buildworker/worker/package_linux64/build/src/interpreter.c:75
eval at /buildworker/worker/package_linux64/build/src/interpreter.c:242
jl_interpret_toplevel_expr at /buildworker/worker/package_linux64/build/src/interpreter.c:34
jl_toplevel_eval_flex at /buildworker/worker/package_linux64/build/src/toplevel.c:577
jl_parse_eval_all at /buildworker/worker/package_linux64/build/src/ast.c:873
jl_load at /buildworker/worker/package_linux64/build/src/toplevel.c:616
include_from_node1 at ./loading.jl:576
unknown function (ip: 0x7f04cb29138b)
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
include at ./sysimg.jl:14
unknown function (ip: 0x7f04cb1299eb)
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
process_options at ./client.jl:305
_start at ./client.jl:371
unknown function (ip: 0x7f04cb29fd28)
jl_call_fptr_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:339 [inlined]
jl_call_method_internal at /buildworker/worker/package_linux64/build/src/julia_internal.h:358 [inlined]
jl_apply_generic at /buildworker/worker/package_linux64/build/src/gf.c:1926
jl_apply at /buildworker/worker/package_linux64/build/ui/../src/julia.h:1424 [inlined]
true_main at /buildworker/worker/package_linux64/build/ui/repl.c:127
main at /buildworker/worker/package_linux64/build/ui/repl.c:264
__libc_start_main at /lib/x86_64-linux-gnu/libc.so.6 (unknown line)
unknown function (ip: 0x4016bc)
unknown function (ip: 0xffffffffffffffff)
Allocations: 527200645 (Pool: 527197857; Big: 2788); GC: 223016
