# Setting Up environment

## Automatic install
Just execute `./runthis.sh`.
Type `sh runthis.sh` inside parent directory.

## Manual Install

### Setup PATH
Edit the `.zshrc` and add `PATH` and `PYTHONPATH` there.

```
export PATH="${HOME}/bin:$PATH"
export PYTHONPATH="${HOME}:pythonlib:${PYTHONPATH}"
```

### Creating directories and copy files
And then create both the directories
```
mkdir $HOME/bin
mkdir $HOME/pythonlib
```

Now copy content of `lib` directory to `pythonlib` in your `HOME` directory.

Copy content of `bin` directory to `bin` in your `HOME` directory.
