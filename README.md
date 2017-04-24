# A Comparison of Genetic Regulatory Network Dynamics and Encoding

This is the corresponding repository for the publication "A Comparison of
Genetic Regulatory Network Dynamics and Encoding", to be presented at the 2017
Genetic and Evolutionary Computation Conference (GECCO).

Video of a GRN performing the Ship Escape problem can be
found [here](https://vimeo.com/214141848).

This work uses both the [GAGA](https://github.com/jdisset/gaga/)
and [grgen](https://github.com/jdisset/grgen) libraries. They were modified for
this work and are included in this repository in their modified forms.

To compile the source code for this experiment, use CMake:

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

A `bin` folder will then be created with `evo`, `flappyKeyboard`, and `dbg`
binaries. `evo` runs GAGA to evolve a GRN and can be launched as follows:

```bash
bin/evo
```

`flappyKeyboard` allows the user to play the FlappyBird game in a virtual
console. The Z key is used to make the bird jump.

```bash
bin/flappyKeyboard
```

The `dbg` binary allows for the execution of a single DNA file, a JSON
representation of a GRN:

```bash
bin/dbg flappy.dna
```

Two example GRNs, `flappy.dna` and `ship.dna` have been included for the Flappy
Bird and Ship Escape problems, respectively.

The Ship Escape problem has been included in headless form but a viewer can be
found at the [ShipEscape](https://github.com/jdisset/shipEscape) repository.
Note that the definition SHIP_DISPLAY must be given in `ship.hpp` to activate
the necessary viewer functions.
