# `TS_find.py`

This is the **only one** necessary for the algorithm to work. It also requires xTB and/or ORCA to be installed. These programs are used to compute the gradient file and for constraint optimisation.

## How to use:
### from terminal

1. Download the project:
```
git clone https://github.com/ER-Elivian/ts_search.git
```
2. Save your structure in the file to_opt.xyz

3. Write bonds_to_search file (see ./tests for examples)

4. Run TS search (need to modify tests/da_test, or leave for test calculation):

```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -tr 100 -p xtb --steps 300 --print
```
where: 
* `-tf` is the force threshold (the maximum force along any bond from bonds_to_search must be less than this value).
* `-tr` is the threshold for the relative excess of the forces on the background (rel excess must be less than this value).
* `-p`, `--programm` is the used software ("orca" or "xtb")
* `--steps`, `-s` is maximum number of steps. the search is interrupted if this value is reached
* `--print` is flag to print output

You can also use ORCA for TS search:
```
python TS_find.py tests/da_test/to_opt.xyz -tf 0.0001 -p orca -s 300 -onp 8 -omm 1500 -OPATH /your/path/to/orca -oms "B3LYP ma-def2-TZVP"
```
<b>NOTE:</b> you may use one or both thresholds

where: 
* `-onp`, `--orca-number-processors` is number of procesors using to calculation
* `-omm`, `--orca-memory` is amount of memory per processor  
* `-OPATH`, `--ORCA-PATH` is ORCA `PATH`. It is necessary if `-onp` value is >1
* `-oms`, `--orca-method-string` is string containing method and basis, it will be written in ORCA's input file after `!` sign

### in-code Python using:

from TS_find import optTS
```
optTS(xyz_path="to_opt.xyz",
      threshold_rel=8,
      threshold_force=0.00004, 
      mirror_coef=1, 
      print_output=True, 
      maxstep=10**4, 
      programm=dict(name="xtb", force_constant= 6))
```
where:
- `xyz_path` is path to xyz file 
- `threshold_rel=8` is threshold for relative excess of forces on  over the background (rel excess must be less this value)
- `threshold_force=0.00004` is threshold for force (maximum force along any bond from bonds_to_search must be less this value)
- `maxstep=10**4` is maximum number of steps. the search is interrupted if this value is reached
- `mirror coef` is the value by which the force projection is multiplied when the longitudinal component relative to the phase vector is reflected. A decrease entails a decrease in velocity, an increase may be the cause of oscillations near the transition state.
- `print_output` is flag to print output

`threshold_rel` or `threshold_force` must exceed 0. Recommended `threshold_rel`>5, `threshold_force`>0.00002, the less, the more accurately the TS will be found, but at the same time, the longer it will take to find it. If botf `threshold_rel` and `threshold_force` is non-zero both thresholds must converged

## There will be detailed description of the algorithm

In general, the convergence of the algorithm on structures with more than 2 bonds has not been proven, but so far no exceptions have been found (exclude cases that have long-periodic oscillation casued by forming unexpeted bonds while being optimized)

. . .


# `test_opt_xtb.py`
checking the convergence of the algorithm in xtb. Prodces table contained RMSD to TS calculated with small constrains (`threshold_force`=5e-7) and number of steps