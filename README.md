# toppic_claude

## Git LFS is required

Some of the runtime resources under `resources/` are large binary/data blobs
and are stored with [Git LFS](https://git-lfs.com/) rather than in the normal
Git history:

- `resources/envcnn_models/*.onnx`, `resources/ecscore_models/*.onnx` — the
  EnvCNN / ECScore neural-network models.
- `resources/base_data/theo_patt.txt` — the theoretical isotope-pattern table.

You **must** have Git LFS installed to get the real files. Without it, a plain
`git clone` leaves small text *pointer* files in their place, and the build /
the tools will fail to load the models and tables.

### Install Git LFS (one time per machine)

```sh
# Debian/Ubuntu
sudo apt-get install git-lfs
# macOS (Homebrew)
brew install git-lfs
# Then register the Git hooks/filters for your user
git lfs install
```

### Clone

With Git LFS installed, a normal clone fetches the LFS files automatically:

```sh
git clone https://github.com/liuxiaowen/toppic_claude.git
```

### Already cloned without Git LFS?

If you cloned before installing Git LFS (so the `.onnx` / `theo_patt.txt`
files are pointer text), install it as above and then pull the real blobs:

```sh
git lfs install
git lfs pull
```
