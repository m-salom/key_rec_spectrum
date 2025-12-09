# Key reconstruction from distance spectrum
This repository provides a C implementation of key reconstruction algorithms using the distance spectrum.

## Prog

### Compilation

Compile without distance multiplicity

    make all

Compile with distance multiplicity

    make all FULL=1

### Usage

    ./keyrec r d n s mode percent_dist_out percent_dist_in

**Arguments**

| Argument | Description |
|---------|-------------|
| `r` | Block length of the code |
| `d` | Block weight of the code |
| `n` | Number of samples |
| `s` | Seed |
| `mode` | Reconstruction mode:<br>• `full` — key reconstruction from a full spectrum<br>• `partial` — reconstruction from a partial spectrum<br>• `both` |
| `percent_dist_out` | % of known distances **outside** the spectrum for partial reconstruction |
| `percent_dist_in` |  % of known distances **inside** the spectrum  for partial reconstruction |

### Example

    ./keyrec 12323 71 10000 7899876 both 70 1
    
---

## Synd Attack

### Create an instance

By default, instance is created in `bike134-12323/` :

    ./build_instance

### Compile the instance

    ./compile bike134-12323

### Run the attack on a random key

### Options

- `-l` : print to ad-hoc file (default: stdout)
- `-v` : verbose level (`-v 7` prints everything)
- `-N` : number of samples (default: 1 000 000)
- `-n` : sigma for noisy syndrome weight experiments
- `-m` : reconstruction mode :
  - `partial` — key reconstruction from a partial spectrum
  - `soft` — key reconstruction from a soft spectrum
  - `both`

### Example

    ./bike134-12323/synd_attack -v 3 -l -N 280000 -n 30 -m soft
