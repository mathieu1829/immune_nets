# immune_nets

# Prerequisites
- [nix package manager](https://nixos.org/download)

Installation

1. run nix-shell setup command:

```bash
nix-shell
```

postgresql will start itself after running this command on the port defined by the `PORT` variable in the `shell.nix` file. All required packages are also downloaded, such as python and sratoolkit are provided by nix-shell and only available in it

2. to download data sets from cellranger execute the command :

```bash
dcrds
```