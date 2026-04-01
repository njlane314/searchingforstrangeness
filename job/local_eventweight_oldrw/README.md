This directory vendors the old-RW NuMI EventWeight FHiCL chain used by the
staged campaign XMLs and local `dev/` validation wrappers.

These files were copied from the corresponding `ubsim/EventWeight/jobs/*`
definitions because the `v08_00_00_82` gpvm stack used for this repo does not
ship the `*_oldflux_rw.fcl` chain in its installed FHiCL path.

Only the missing EventWeight configuration chain is mirrored here. Standard
MicroBooNE service includes such as `services_microboone.fcl` still resolve
from the normal `uboonecode` environment.
