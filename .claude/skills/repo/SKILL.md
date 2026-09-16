---
name: repo
description: Publish the private TopPIC Suite repository (liuxiaowen/toppic_suite_private) to the public GitHub repository (toppic-suite/toppic-suite) - sync the main branch and the release tags, carry the Git LFS model files across, and handle the public branch's diverged history. Use when asked to update, sync, publish or push the public repo, or to make a tagged release public.
---

# Updating the public repository from the private one

Two GitHub repositories hold the TopPIC Suite:

| Role | Repository | Branch |
|---|---|---|
| **Private, the source of truth** | `https://github.com/liuxiaowen/toppic_suite_private.git` (the remote of this checkout may still show its old name `proteomics_cpp`; GitHub redirects) | `main` |
| **Public release copy** | `https://github.com/toppic-suite/toppic-suite.git` | `main` |

Development happens in the private repository. Publishing means making the
public `main` carry the private `main`'s content plus all tags. Both default
branches are called `main`; neither repository has a `master`.

The public repository is **not** a plain mirror: it has merged at least one
external pull request (`ece41d25`, "Supported Mac build") that the private
history lacks, so the two `main` branches have diverged. A plain
`git merge` between them conflicts (CMakeLists.txt, README.md,
file_util.cpp, and a file the private tree deleted). The procedure below
keeps the public history but takes the private tree verbatim, which is what
a sync should produce, and works without manual conflict resolution.

## Prerequisites

- `git` and `git-lfs`, with `git lfs install` run once on the machine. Three
  files under `res/` (two `.onnx` models, `theo_patt.txt`, about 105 MB) are
  LFS objects. **LFS objects are stored per repository**: a `git fetch` only
  brings their pointers, so they must be fetched from the private
  repository and pushed to the public one explicitly (steps below).
- Push rights on `toppic-suite/toppic-suite` (an HTTPS token or SSH key that
  git can use non-interactively).
- The private `main` is at the state to publish (merge `release_1.9` into
  `main` first if needed) and the release tag exists there.

## Procedure

Run in a throwaway directory; nothing here touches this checkout.

```sh
set -e
work=$(mktemp -d)
cd "$work"
git clone https://github.com/toppic-suite/toppic-suite.git toppic_suite
cd toppic_suite
git remote add upstream https://github.com/liuxiaowen/toppic_suite_private.git
git fetch upstream --tags
git lfs fetch upstream main            # the model files behind the pointers

# 1. Does the public main have commits the private one lacks?
git log --oneline upstream/main..main
```

If that list is **empty**, fast-forward:

```sh
git merge --ff-only upstream/main
```

If it is **not empty** (the normal case today), first decide whether any of
those commits carry changes that should survive; if so, port them into the
private repository (cherry-pick, then push there and start over). Then
merge keeping the public history but taking the private tree exactly:

```sh
git merge -s ours --no-commit upstream/main   # merge commit, both parents
git read-tree -u --reset upstream/main        # ...with the private tree
git commit -m "Merge private main into public main"
git diff --stat upstream/main HEAD            # must print nothing
```

Then publish:

```sh
git lfs push origin main                      # LFS objects first
git push origin main
git push origin --tags
cd /
rm -rf "$work"
```

Verify from anywhere:

```sh
git ls-remote --heads https://github.com/toppic-suite/toppic-suite.git main
git ls-remote --tags  https://github.com/toppic-suite/toppic-suite.git v1.9.0.0
```

The public `main` hash is the new merge commit (or, after a fast-forward,
the private `main` hash), and the tag resolves to the same commit as in the
private repository.

## Optional: publish other branches

Only `main` and the tags are synced above. The public repository also lacks
`release_1.9` (and the private working branches `mdec`, `post_mass_match`,
`low_resolution`, which stay private). To publish the release branch too:

```sh
git push origin upstream/release_1.9:refs/heads/release_1.9
```

(from inside the throwaway clone, after `git fetch upstream`). LFS objects
on that branch are the same as on `main` once it has been merged there.

## Why the original one-liner script was changed

The first version of this procedure was:

```sh
git clone .../toppic-suite.git toppic_suite && cd toppic_suite
git remote add upstream .../toppic_suite_private.git
git pull; git fetch upstream --tags; git merge upstream/master
git push --tags; git push
```

It fails for four reasons, now fixed above:

1. `upstream/master` does not exist; the branch is `main`.
2. `git merge upstream/main` stops with conflicts in four files, after
   which the pushes push nothing new and the cleanup deletes the half-merged
   clone. The `-s ours` + `read-tree` merge never conflicts.
3. Nothing fetched or pushed the LFS objects, so the public checkout would
   have contained pointer files whose objects do not exist on that server
   ("Object does not exist on the server" at clone or checkout time).
4. `git pull` right after `git clone` is a no-op, and without `set -e` a
   failure in one step did not stop the following ones.

## Troubleshooting

- **`! [rejected] ... (would clobber existing tag)`** on `git fetch
  upstream --tags`: a tag was moved in the private repository (this was
  done for `v1.9.0.0` once). If the private tag is the intended one, run
  `git fetch upstream --tags --force` and later `git push origin --tags
  --force`; announce it, since existing clones keep the old tag.
- **`Object does not exist on the server`** or `LFS: object not found`
  during `git push`: an LFS object referenced by a pushed commit is not in
  the local cache. `git lfs fetch upstream --all` (every version of every
  LFS file; several hundred MB) then `git lfs push origin --all`.
- **`git diff --stat upstream/main HEAD` prints files** after the merge:
  the read-tree step was skipped or ran against the wrong ref; do not
  push. `git reset --hard main` (the clone's original tip) and redo it.
- **Authentication prompts in a script**: configure a credential helper or
  use an SSH URL for `origin`; the private repository is read-only here and
  needs no push rights.
- **Never force-push `main`** to make the public repository a mirror: it
  would erase the contributor's merged pull request from the public
  history. The merge above keeps it as an ancestor while the tree follows
  the private repository.
